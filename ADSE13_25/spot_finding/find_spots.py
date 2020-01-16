from __future__ import absolute_import, division, print_function

import logging
import os, glob
logger = logging.getLogger('exafel.find_spots')
from dxtbx.model.experiment_list import ExperimentListFactory, ExperimentList, ExperimentListDumper
from dials.algorithms.spot_prediction import StillsReflectionPredictor
from dials.array_family import flex
from libtbx.utils import Abort, Sorry
from libtbx.phil import parse


# Imports for LS49
from dials.command_line.stills_process import Script, Processor, control_phil_str, dials_phil_str, program_defaults_phil_str

LS49_phil_str = """
LS49 {
  dump_CBF = False
    .type = bool
    .help = if True, then dumps images as CBFs. Can specificy a list of timestamps to dump
  predict_spots = False
    .type = bool
    .help = If False, only do spotfinding. Ideally used for Jungfrau XTC data. If True, also does downstream analysis \
            like spot prediction using crystal model indexed on rayonix.
  path_to_rayonix_crystal_models = None
    .type = str
    .help = path to crystal models that were indexed using the rayonix data.Should only be a path \
            Eg. /my/path/to/rayonix_models/ \
            Where the folder contains all the integrated_experiments.json files from the rayonix analysis \
            Note that this option should be provided if predict_spots is True \ 
            Update 15th Jan, 2020 --> Also using this to figure out timestamps for dumping CBFs \
            Folder should contain list of int-x files
  path_to_jungfrau_detector_model = None
    .type = str
    .help = path to jungfrau detector models that were indexed. Should be path to the experiment file including the file itself\
            Eg. /my/path/to/jungfrau_model/experiment.json \
            Note that this option should be provided if predict_spots is True
}
"""
phil_scope = parse(control_phil_str + dials_phil_str + LS49_phil_str, process_includes=True).fetch(parse(program_defaults_phil_str))


message = ''' This is a specialized version of stills_process which only does spot_finding. Difference is it only
              prints out the datablock/strong pickle files of those images which have hits on them
'''


def do_import(filename, load_models=True):
    logger.info("Loading %s" % os.path.basename(filename))
    experiments = ExperimentListFactory.from_filenames([filename], load_models=False)
    if len(experiments) == 0:
        try:
            experiments = ExperimentListFactory.from_json_file(filename)
        except ValueError:
            raise Abort("Could not load %s" % filename)

    if len(experiments) == 0:
        raise Abort("Could not load %s" % filename)

    from dxtbx.imageset import ImageSetFactory

    for experiment in experiments:
        if load_models:
            experiment.load_models()
        imageset = ImageSetFactory.imageset_from_anyset(experiment.imageset)
        imageset.set_scan(None)
        imageset.set_goniometer(None)
        experiment.imageset = imageset
        experiment.scan = None
        experiment.goniometer = None

    return experiments

class SpotFinding_Script(Script):
    def __init__(self):
        """Initialise the script."""
        from dials.util.options import OptionParser
        import libtbx.load_env

        # The script usage
        usage = (
            "usage: %s [options] [param.phil] filenames" % libtbx.env.dispatcher_name
        )

        self.tag = None
        self.reference_detector = None

        # Create the parser
        self.parser = OptionParser(usage=usage, phil=phil_scope, epilog="")

    def run(self):
        """Execute the script."""
        from dials.util import log
        from time import time
        from libtbx import easy_mp
        import copy

        # Parse the command line
        params, options, all_paths = self.parser.parse_args(
            show_diff_phil=False, return_unhandled=True, quick_parse=True
        )

        # Check we have some filenames
        if not all_paths:
            self.parser.print_help()
            return

        # Mask validation
        for mask_path in params.spotfinder.lookup.mask, params.integration.lookup.mask:
            if mask_path is not None and not os.path.isfile(mask_path):
                raise Sorry("Mask %s not found" % mask_path)

        # Save the options
        self.options = options
        self.params = params

        st = time()

        # Configure logging
        #log.config(
        #    params.verbosity, info="exafel_spotfinding.process.log", debug="exafel.spot_finding.debug.log"
        #)

        bad_phils = [f for f in all_paths if os.path.splitext(f)[1] == ".phil"]
        if len(bad_phils) > 0:
            self.parser.print_help()
            logger.error(
                "Error: the following phil files were not understood: %s"
                % (", ".join(bad_phils))
            )
            return

        # Log the diff phil
        diff_phil = self.parser.diff_phil.as_str()
        if diff_phil is not "":
            logger.info("The following parameters have been modified:\n")
            logger.info(diff_phil)

        for abs_params in self.params.integration.absorption_correction:
            if abs_params.apply:
                if not (
                    self.params.integration.debug.output
                    and not self.params.integration.debug.separate_files
                ):
                    raise Sorry(
                        "Shoeboxes must be saved to integration intermediates to apply an absorption correction. "
                        + "Set integration.debug.output=True, integration.debug.separate_files=False and "
                        + "integration.debug.delete_shoeboxes=True to temporarily store shoeboxes."
                    )

        self.load_reference_geometry()
        from dials.command_line.dials_import import ManualGeometryUpdater

        update_geometry = ManualGeometryUpdater(params)

        # Import stuff
        logger.info("Loading files...")
        pre_import = params.dispatch.pre_import or len(all_paths) == 1
        if True: #pre_import:
            # Handle still imagesets by breaking them apart into multiple experiments
            # Further handle single file still imagesets (like HDF5) by tagging each
            # frame using its index
            experiments = ExperimentList()
            for path in all_paths:
                experiments.extend(do_import(path, load_models=False))

            indices = []
            basenames = []
            split_experiments = []
            for i, imageset in enumerate(experiments.imagesets()):
                assert len(imageset) == 1
                paths = imageset.paths()
                indices.append(i)
                basenames.append(os.path.splitext(os.path.basename(paths[0]))[0])
                split_experiments.append(experiments[i : i + 1])
            tags = []
            for i, basename in zip(indices, basenames):
                if basenames.count(basename) > 1:
                    tags.append("%s_%05d" % (basename, i))
                else:
                    tags.append(basename)

            # Wrapper function
            def do_work(i, item_list):
                processor = SpotFinding_Processor(copy.deepcopy(params), composite_tag = "%04d"%i, rank = i)
                if params.LS49.dump_CBF:
                    print ('READING IN TIMESTAMPS TO DUMP')
                    # Read in file with timestamps information
                    processor.timestamps_to_dump = []
                    for fin in glob.glob(os.path.join(self.params.LS49.path_to_rayonix_crystal_models, 'int-0-*')):
                      int_file=os.path.basename(fin)
                      ts = int_file[6:23] 
                      processor.timestamps_to_dump.append(ts)
                    #with open(os.path.join(self.params.output.output_dir,'../timestamps_to_dump.dat'), 'r') as fin:
                    #    for line in fin:
                    #        if line !='\n':
                    #            ts = line.split()[0].strip()
                    #            processor.timestamps_to_dump.append(ts)

                from dials.array_family import flex
                all_spots_from_rank = flex.reflection_table()
                for item in item_list:
                    try:
                        assert len(item[1]) == 1
                        experiment = item[1][0]
                        experiment.load_models()
                        imageset = experiment.imageset
                        update_geometry(imageset)
                        experiment.beam = imageset.get_beam()
                        experiment.detector = imageset.get_detector()
                    except RuntimeError as e:
                        logger.warning(
                          "Error updating geometry on item %s, %s"
                          % (str(item[0]), str(e)))
                        continue

                    if self.reference_detector is not None:
                        from dxtbx.model import Detector
                        experiment = item[1][0]
                        imageset = experiment.imageset
                        imageset.set_detector(Detector.from_dict(self.reference_detector.to_dict()))
                        experiment.detector = imageset.get_detector()

                    refl_table=processor.process_experiments(item[0], item[1], item[2])
                    if refl_table is not None:
                        all_spots_from_rank.extend(refl_table)
                processor.finalize()
                return all_spots_from_rank

            iterable = zip(tags, split_experiments, indices)

        # Process the data
        if True: #params.mp.method == 'mpi':
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            rank = comm.Get_rank() # each process in MPI has a unique id, 0-indexed
            size = comm.Get_size() # size: number of processes running in this job

        # Configure the logging
            if params.output.logging_dir is None:
                info_path = ''
                debug_path = ''
            else:
                import sys
                log_path = os.path.join(params.output.logging_dir, "log_rank%04d.out"%rank)
                error_path = os.path.join(params.output.logging_dir, "error_rank%04d.out"%rank)
                print ("Redirecting stdout to %s"%log_path)
                print ("Redirecting stderr to %s"%error_path)
                sys.stdout = open(log_path,'a', buffering=0)
                sys.stderr = open(error_path,'a',buffering=0)
                print ("Should be redirected now")

                info_path = os.path.join(params.output.logging_dir, "info_rank%04d.out"%rank)
                debug_path = os.path.join(params.output.logging_dir, "debug_rank%04d.out"%rank)

            from dials.util import log
            print ('IOTA_ALL_SPOTS_RANKS_0')
            #log.config(params.verbosity, info=info_path, debug=debug_path)
            subset = [item for i, item in enumerate(iterable) if (i+rank)%size == 0]
            all_spots_from_rank = do_work(rank, subset)
            all_spots_rank0 = comm.gather(all_spots_from_rank, root=0)
            print ('IOTA_ALL_SPOTS_RANKS_1')
            exit()
            if rank == 0:
                from dials.array_family import flex
                all_spots = flex.reflection_table()
                for ii,refl_table in enumerate(all_spots_rank0):
                    if refl_table is not None:
                        all_spots.extend(refl_table)
                from libtbx.easy_pickle import dump
                #dump('all_spots.pickle', all_spots_rank0)
                #dump('all_experiments.pickle', experiments)
                #print ('IOTA_ALL_SPOTS_RANKS_2')
                #print ('IOTA_ALL_SPOTS_RANKS_3')
                from dials.algorithms.spot_finding import per_image_analysis
                from six.moves import cStringIO as StringIO
                s = StringIO()
                # Assuming one datablock. Might be dangerous
                # FIXME
                from dxtbx.format.cbf_writer import FullCBFWriter
                for i, imageset in enumerate(experiments.imagesets()):
                    print("Number of centroids per image for imageset %i:" %i, file=s)
                    #from IPython import embed; embed(); exit()
                    print ('IOTA_ALL_SPOTS_RANKS_4')
                    stats = custom_stats_imageset(
                        imageset, all_spots.select(all_spots['img_id'] == i))
                    n_spots_total = flex.int(stats.n_spots_total)
                    max_number_of_spots = max(stats.n_spots_total)
                    for num_spots in range(1,max_number_of_spots+1):
                        print ("IOTA_NUMBER_OF_SPOTS %d %d"%(num_spots, len(n_spots_total.select(n_spots_total==num_spots))))
                    if max_number_of_spots > 0:
                        # assuming one imageset per experiment here : applicable for stills
                        ts = imageset.get_image_identifier(0)
                        xfel_ts = ts[0:4] + ts[5:7] + ts[8:10] + ts[11:13] + ts[14:16] + ts[17:19] + ts[20:23] 
                        cbf_path = os.path.join(params.output.logging_dir, 'jungfrau_%s.cbf'%xfel_ts)
                        cbf_writer = FullCBFWriter(imageset=imageset)
                        cbf_writer.write_cbf(cbf_path)
                    per_image_analysis.print_table(stats)
                    logger.info(s.getvalue())
            comm.barrier() 

class SpotFinding_Processor(Processor):
    def process_experiments(self, tag, experiments, img_id):
        import os
        if not self.params.output.composite_output:
            self.setup_filenames(tag)
        self.img_id = img_id
        self.tag = tag
        self.debug_start(tag)

#    if not self.params.output.composite_output and self.params.output.datablock_filename:
#      from dxtbx.datablock import DataBlockDumper
#      dump = DataBlockDumper(datablock)
#      dump.as_json(self.params.output.datablock_filename)


        try:
            if self.params.LS49.dump_CBF:
                from dxtbx.format.cbf_writer import FullCBFWriter
            # assuming one imageset per experiment here : applicable for stills
                ts = experiments[0].imageset.get_image_identifier(0)
                xfel_ts = ts[0:4] + ts[5:7] + ts[8:10] + ts[11:13] + ts[14:16] + ts[17:19] + ts[20:23]
                print ('TIMESTAMPS = ',xfel_ts)
                if xfel_ts in self.timestamps_to_dump:
                    cbf_path = os.path.join(self.params.output.logging_dir, 'jungfrauhit_%s.cbf'%xfel_ts)
                    cbf_writer = FullCBFWriter(imageset=experiments[0].imageset)
                    cbf_writer.write_cbf(cbf_path)
                return None
        except Exception as e:
            print ('Error dumping CBFs', tag, str(e))
    # Do spotfinding
        try:
            self.debug_write("spotfind_start")
            observed = self.find_spots(experiments)
            if not self.params.LS49.predict_spots:
                return observed
        except Exception as e:
            print("Error spotfinding", tag, str(e))
            return None


        try:
            self.debug_write('spot_prediction_start')
            observed = self.predict_spots_from_rayonix_crystal_model(experiments, observed)
            return observed
        except Exception as e:
            print("Error spotfinding - in spot_prediction", tag, str(e))
            return None

    def predict_spots_from_rayonix_crystal_model(self, experiments, observed):
        """ Reads in the indexed rayonix model, predict spots using the crystal model on the jungfrau detector"""
        pass
        # Make sure experimental model for rayonix is supplied. Also the experimental geometry of the jungfrau is supplied
        assert self.params.LS49.path_to_rayonix_crystal_models is not None, 'Rayonix crystal model path is empty. Needs to be specified'
        assert self.params.LS49.path_to_jungfrau_detector_model is not None, 'Jungfrau_detector model path is empty. Needs to be specified'
        ts = self.tag.split('_')[-1] # Assuming jungfrau cbfs are names as 'jungfrauhit_20180501133315870' 
        # Load rayonix experimental model
        rayonix_fname=os.path.join(self.params.LS49.path_to_rayonix_crystal_models, 'idx-%s_integrated_experiments.json'%ts)
        rayonix_expt=ExperimentListFactory.from_json_file(rayonix_fname, check_format=False)
        jungfrau_det=ExperimentListFactory.from_json_file(self.params.LS49.path_to_jungfrau_detector_model, check_format=False)
        # Reset stuff here
        # Should have 
        # a. Jungfrau detector geometry
        # b. Rayonix indexed crystal model
        from dials.algorithms.refinement.prediction.managed_predictors import ExperimentsPredictorFactory
        from dials.algorithms.indexing import index_reflections
        experiments[0].detector=jungfrau_det[0].detector
        experiments[0].crystal=rayonix_expt[0].crystal
        if False:
            observed['id']=flex.int(len(observed), -1)
            observed['imageset_id']=flex.int(len(observed), 0)
            observed.centroid_px_to_mm(experiments[0].detector, experiments[0].scan)
            observed.map_centroids_to_reciprocal_space(experiments[0].detector, experiments[0].beam, experiments[0].goniometer)
            index_reflections(observed, experiments)
            ref_predictor=ExperimentsPredictorFactory.from_experiments(experiments)
            ref_predictor(observed)
            observed['id']=flex.int(len(observed), 0)
            from libtbx.easy_pickle import dump
            dump('my_observed_prediction_%s.pickle'%self.tag, observed)
            dumper=ExperimentListDumper(experiments)
            dumper.as_json('my_observed_prediction_%s.json'%self.tag)
        
        predictor=StillsReflectionPredictor(experiments[0])
        ubx=predictor.for_ub(experiments[0].crystal.get_A())
        ubx['id']=flex.int(len(ubx),0)
        n_predictions = len(ubx)
        n_observed = len(observed)
        if len(observed) > 3 and len(ubx) >= len(observed):
            from libtbx.easy_pickle import dump
            dump('my_prediction_%s.pickle'%self.tag, ubx)
            dumper=ExperimentListDumper(experiments)
            dumper.as_json('my_prediction_%s.json'%self.tag)
            #from IPython import embed; embed(); exit()
            exit()
        
        

    def find_spots(self, experiments):
        from time import time
        from dials.array_family import flex
        st = time()

        logger.info('*' * 80)
        logger.info('Finding Strong Spots')
        logger.info('*' * 80)

        # Find the strong spots
        observed = flex.reflection_table.from_observations(experiments, self.params)

        # Reset z coordinates for dials.image_viewer; see Issues #226 for details
        xyzobs = observed['xyzobs.px.value']
        for i in xrange(len(xyzobs)):
            xyzobs[i] = (xyzobs[i][0], xyzobs[i][1], 0)
            bbox = observed['bbox']
        for i in xrange(len(bbox)):
            bbox[i] = (bbox[i][0], bbox[i][1], bbox[i][2], bbox[i][3], 0, 1)

        if self.params.output.composite_output:
            pass # no composite strong pickles yet
        else:
            # Save the reflections to file
            # Only save those which have spots
            logger.info('\n' + '-' * 80)
            if self.params.output.strong_filename and len(observed) > 0:
                self.save_reflections(observed, self.params.output.strong_filename)
                from dxtbx.model.experiment_list import ExperimentListDumper
                dump = ExperimentListDumper(experiments)
                dump.as_json(self.params.output.experiments_filename)

        logger.info('')
        logger.info('Time Taken = %f seconds' % (time() - st))
        observed['img_id'] = flex.int(len(observed), self.img_id)
        return observed


def custom_stats_imageset(imageset, reflections, resolution_analysis=False, plot=False):
    from dials.algorithms.spot_finding import per_image_analysis
    from libtbx import group_args
    from dials.array_family import flex
    n_spots_total = []
    n_spots_no_ice = []
    n_spots_4A = []
    total_intensity = []
    estimated_d_min = []
    d_min_distl_method_1 = []
    d_min_distl_method_2 = []
    noisiness_method_1 = []
    noisiness_method_2 = []

    try:
        start, end = imageset.get_array_range()
    except AttributeError:
        start = 0
    for i in range(len(imageset)):
        stats = per_image_analysis.stats_single_image(
          imageset[i:i+1],
          reflections, i=i+start,
          #reflections.select(reflections['img_id']==i+start), i=i+start,
          resolution_analysis=resolution_analysis, plot=plot)
        n_spots_total.append(stats.n_spots_total)
        n_spots_no_ice.append(stats.n_spots_no_ice)
        n_spots_4A.append(stats.n_spots_4A)
        total_intensity.append(stats.total_intensity)
        estimated_d_min.append(stats.estimated_d_min)
        d_min_distl_method_1.append(stats.d_min_distl_method_1)
        noisiness_method_1.append(stats.noisiness_method_1)
        d_min_distl_method_2.append(stats.d_min_distl_method_2)
        noisiness_method_2.append(stats.noisiness_method_2)

    return group_args(n_spots_total=n_spots_total,
                      n_spots_no_ice=n_spots_no_ice,
                      n_spots_4A=n_spots_4A,
                      total_intensity=total_intensity,
                      estimated_d_min=estimated_d_min,
                      d_min_distl_method_1=d_min_distl_method_1,
                      noisiness_method_1=noisiness_method_1,
                      d_min_distl_method_2=d_min_distl_method_2,
                      noisiness_method_2=noisiness_method_2)




if __name__ == '__main__':
    #with dials.util.show_mail_on_error():
    try:
        script = SpotFinding_Script()
        script.run()
    except Exception as e:
        print (e)
