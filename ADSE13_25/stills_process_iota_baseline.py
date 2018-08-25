from __future__ import absolute_import, division, print_function

import logging
import os

from dxtbx.datablock import DataBlockFactory
from libtbx.utils import Abort, Sorry
from libtbx.phil import parse

from dials.command_line.stills_process import control_phil_str, dials_phil_str, program_defaults_phil_str
from dials.command_line.stills_process import do_import, Script, Processor

logger = logging.getLogger('stills_process_iota')

help_message = '''
script for processing stills using IOTA-style to explore optimal spotfinding/indexing params
'''

iota_phil_str = '''
iota {
  method = off *random_sub_sampling
    .type = choice
    .help = Type of IOTA processing to be done. \
            off : No IOTA processing is done. \
            random-sub-sampling : randomly sub-sample observed bragg spots and index. Can be done multiple times. See options for random-sub-sampling if this is used.
  random_sub_sampling {
    ntrials = 10
      .type = int
      .help = Number of random sub-samples to be selected
    fraction_sub_sample = 0.2
      .type = float
      .help = fraction of sample to be sub-sampled. Should be between 0 and 1
    consensus_function = *unit_cell
      .type = choice
      .help = choose the type of consensus function to be employed for random_sub_sampling. More details \
              in the functions themselves
    show_plot = False
      .type = bool
      .help = Flag to indicate whether plots for clustering are to be displayed. Useful for debugging
    no_outlier_rejection_and_candidates_refinement=False
      .type = bool
      .help = Flag to indicate if candidate basis vectors should be refined and whether \
              outlier rejectionis needed
    finalize_method = union_and_reindex *reindex_with_known_crystal_models
      .type = choice
      .help = union_and_reindex will take union of all spots used to obtain \
              cluster (hence lattice) and then reindex with all the lattice models \
              of that cluster.\
              reindex_with_known_crystal_models will just index the spots with \
              known_orientation_indexer. Useful if clustering fails but indexing \
              succeeds in limited trials
    Z_cutoff = 1.0
      .type = float
      .help = Z-score cutoff for accepting/rejecting bragg spots based on difference between \
              fractional and integer hkl. This will be used for finalize_method = union_and_reindex

  }
}
'''

phil_scope = parse(control_phil_str + dials_phil_str + iota_phil_str, process_includes=True).fetch(parse(program_defaults_phil_str))


class Script_iota(Script):
  ''' Script class with functions customized for iota style processsing '''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # The script usage
    usage = "usage: %s [options] [param.phil] filenames" % libtbx.env.dispatcher_name

    self.tag = None
    self.reference_detector = None

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message
      )

  def run(self):
    from dials.util import log
    from time import time
    from libtbx import easy_mp
    import copy

    #Parse command line
    params, options, all_paths = self.parser.parse_args(show_diff_phil=False, return_unhandled=True, quick_parse=True)

    # Check we have some filenames
    if not all_paths:
      self.parser.print_help()
      return

    # Mask validation
    for mask_path in params.spotfinder.lookup.mask, params.integration.lookup.mask:
      if mask_path is not None and not os.path.isfile(mask_path):
        raise Sorry("Mask %s not found"%mask_path)

    # Save the options
    self.options = options
    self.params = params

    st = time()

    # Configure logging
    log.config(
      params.verbosity,
      info='exafel.iota.process.log',
      debug='exafel.iota.process.debug.log')

    # Log the diff phil
    diff_phil = self.parser.diff_phil.as_str()
    if diff_phil is not '':
      logger.info('The following parameters have been modified:\n')
      logger.info(diff_phil)


    for abs_params in self.params.integration.absorption_correction:
      if abs_params.apply:
        if not (self.params.integration.debug.output and not self.params.integration.debug.separate_files):
          raise Sorry('Shoeboxes must be saved to integration intermediates to apply an absorption correction. '\
            +'Set integration.debug.output=True, integration.debug.separate_files=False and '\
            +'integration.debug.delete_shoeboxes=True to temporarily store shoeboxes.')

    self.load_reference_geometry()
    from dials.command_line.dials_import import ManualGeometryUpdater
    update_geometry = ManualGeometryUpdater(params)

    # Import stuff
    logger.info("Loading files...")
    pre_import = params.dispatch.pre_import or len(all_paths) == 1
    if pre_import:
      # Handle still imagesets by breaking them apart into multiple datablocks
      # Further handle single file still imagesets (like HDF5) by tagging each
      # frame using its index

      datablocks = [do_import(path) for path in all_paths]

      indices = []
      basenames = []
      split_datablocks = []
      for datablock in datablocks:
        for imageset in datablock.extract_imagesets():
          paths = imageset.paths()
          for i in xrange(len(imageset)):
            subset = imageset[i:i+1]
            split_datablocks.append(DataBlockFactory.from_imageset(subset)[0])
            indices.append(i)
            basenames.append(os.path.splitext(os.path.basename(paths[i]))[0])
      tags = []
      for i, basename in zip(indices, basenames):
        if basenames.count(basename) > 1:
          tags.append("%s_%05d"%(basename, i))
        else:
          tags.append(basename)

      # Wrapper function
      def do_work(i, item_list):
        print ('Hello pre-import from rank %d'%i)
        processor = Processor_iota(copy.deepcopy(params), composite_tag = "%04d"%i)

        for item in item_list:
          try:
            for imageset in item[1].extract_imagesets():
              update_geometry(imageset)
          except RuntimeError as e:
            logger.warning("Error updating geometry on item %s, %s"%(str(item[0]), str(e)))
            continue

          if self.reference_detector is not None:
            from dxtbx.model import Detector
            for i in range(len(imageset)):
              imageset.set_detector(
                Detector.from_dict(self.reference_detector.to_dict()),
                index=i)

          processor.process_datablock(item[0], item[1])
        processor.finalize()

      iterable = zip(tags, split_datablocks)

    else:
      basenames = [os.path.splitext(os.path.basename(filename))[0] for filename in all_paths]
      tags = []
      for i, basename in enumerate(basenames):
        if basenames.count(basename) > 1:
          tags.append("%s_%05d"%(basename, i))
        else:
          tags.append(basename)

      # Wrapper function
      def do_work(i, item_list):
        processor = Processor_iota(copy.deepcopy(params), composite_tag = "%04d"%i)
        for item in item_list:
          tag, filename = item

          datablock = do_import(filename)
          imagesets = datablock.extract_imagesets()
          if len(imagesets) == 0 or len(imagesets[0]) == 0:
            logger.info("Zero length imageset in file: %s"%filename)
            return
          if len(imagesets) > 1:
            raise Abort("Found more than one imageset in file: %s"%filename)
          if len(imagesets[0]) > 1:
            raise Abort("Found a multi-image file. Run again with pre_import=True")

          try:
            update_geometry(imagesets[0])
          except RuntimeError as e:
            logger.warning("Error updating geometry on item %s, %s"%(tag, str(e)))
            continue

          if self.reference_detector is not None:
            from dxtbx.model import Detector
            imagesets[0].set_detector(Detector.from_dict(self.reference_detector.to_dict()))

          processor.process_datablock(tag, datablock)
        processor.finalize()

      iterable = zip(tags, all_paths)

    # Process the data
    if params.mp.method == 'mpi':
      from mpi4py import MPI
      comm = MPI.COMM_WORLD
      rank = comm.Get_rank() # each process in MPI has a unique id, 0-indexed
      size = comm.Get_size() # size: number of processes running in this job
      subset = [item for i, item in enumerate(iterable) if (i+rank)%size== 0] 
      do_work(rank, subset)
      comm.barrier()
    else:
      from dxtbx.command_line.image_average import splitit
      if params.mp.nproc == 1:
        do_work(0, iterable)
      else:
        result = list(easy_mp.multi_core_run(
          myfunction=do_work,
          argstuples=list(enumerate(splitit(iterable, params.mp.nproc))),
          nproc=params.mp.nproc))
        error_list = [r[2] for r in result]
        if error_list.count(None) != len(error_list):
          print("Some processes failed excecution. Not all images may have processed. Error messages:")
          for error in error_list:
            if error is None: continue
            print(error)

    # Total Time
    logger.info("")
    logger.info("Total Time Taken = %f seconds" % (time() - st))

class Processor_iota(Processor):
  ''' Processor class with functions customized for iota style processing '''
  
  def process_datablock(self, tag, datablock):
    if not self.params.output.composite_output:
      self.setup_filenames(tag)
    self.tag = tag
    print('MP method = ',self.params.mp.method)
    if self.params.output.datablock_filename:
      from dxtbx.datablock import DataBlockDumper
      dump = DataBlockDumper(datablock)
      dump.as_json(self.params.output.datablock_filename)
  
    # Do the processing
    try:
      self.pre_process(datablock)
    except Exception as e:
      print("Error in pre-process", tag, str(e))
      if not self.params.dispatch.squash_errors: raise
      return
    try:
      if self.params.dispatch.find_spots:
        observed = self.find_spots(datablock)
      else:
        print("Spot Finding turned off. Exiting")
        return
    except Exception as e:
      print("Error spotfinding", tag, str(e))
      if not self.params.dispatch.squash_errors: raise
      return
    try:
      if self.params.dispatch.index:
        if self.params.iota.method == 'random_sub_sampling':
          from scitbx.array_family import flex
          len_max_indexed = -999
          experiments_list = []
          # Add an id for each strong spot observed in the image
          observed['spot_id'] = flex.size_t(range(len(observed)))
          # No outlier rejection or refinement should be done for the candidate basis vectors
          outlier_rejection_flag=self.params.indexing.stills.candidate_outlier_rejection
          refine_all_candidates_flag=self.params.indexing.stills.refine_all_candidates
          if self.params.iota.random_sub_sampling.no_outlier_rejection_and_candidates_refinement:
            self.params.indexing.stills.candidate_outlier_rejection=False
            self.params.indexing.stills.refine_all_candidates=False

          observed_samples_list = []
          for trial in range(self.params.iota.random_sub_sampling.ntrials):
            flex.set_random_seed(trial+1001)
            observed_sample = observed.select(flex.random_selection(len(observed), int(len(observed)*self.params.iota.random_sub_sampling.fraction_sub_sample)))
            try:
              print ('IOTA: SUM_INTENSITY_VALUE',sum(observed_sample['intensity.sum.value']), ' ',trial)
              if self.params.iota.random_sub_sampling.finalize_method == 'union_and_reindex':
                experiments_tmp, indexed_tmp = self.index_with_iota(datablock, observed_sample)
              elif self.params.iota.random_sub_sampling.finalize_method == 'reindex_with_known_crystal_models':
                experiments_tmp, indexed_tmp = self.index(datablock, observed_sample)

              experiments_list.append(experiments_tmp)
              observed_samples_list.append(observed_sample)
            except:
              print('Indexing failed for some reason')
          if self.params.iota.random_sub_sampling.consensus_function == 'unit_cell':
            from exafel_project.ADSE13_25.consensus_functions import get_uc_consensus as get_consensus
            known_crystal_models, clustered_experiments_list = get_consensus(experiments_list, show_plot=self.params.iota.random_sub_sampling.show_plot, return_only_first_indexed_model=False, finalize_method=self.params.iota.random_sub_sampling.finalize_method)
          print ('IOTA: Finalizing consensus')
          if self.params.iota.random_sub_sampling.finalize_method == 'reindex_with_known_crystal_models':
            print ('IOTA: Chosen finalize method is reindex_with_known_crystal_models')
            self.known_crystal_models = known_crystal_models
            # Set back whatever PHIL parameter was supplied by user for outlier rejection and refinement
            self.params.indexing.stills.candidate_outlier_rejection=outlier_rejection_flag
            self.params.indexing.stills.refine_all_candidates=refine_all_candidates_flag
            experiments, indexed = self.index(datablock, observed)
            print('fraction subsampled = %5.2f with %d indexed spots ' %(self.params.iota.random_sub_sampling.fraction_sub_sample,len(indexed)))

          elif self.params.iota.random_sub_sampling.finalize_method == 'union_and_reindex':
            print ('IOTA: Chosen finalize method is union_and_reindex')
            # Take union of all spots used to index each lattice cluster
            from dials.array_family import flex as dials_flex
            from dxtbx.model.experiment_list import ExperimentList, Experiment
            indexed = dials_flex.reflection_table()
            experiments = ExperimentList()
            sample = {} 
            all_experimental_models = {}
            for idx,crystal_model in enumerate(clustered_experiments_list):
              if crystal_model >= 0:
                if crystal_model not in sample:
                  sample[crystal_model] = []
                  all_experimental_models[crystal_model] = []
                sample[crystal_model].append(observed_samples_list[idx]['spot_id'])
                all_experimental_models[crystal_model].append(experiments_list[idx])
            for crystal_model in sample: 
              # Need to have a minimum number of experiments for correct stats
              # FIXME number should not be hardcoded. ideally a phil param
              if len(all_experimental_models[crystal_model]) < 3: 
                continue
              self.known_crystal_models = None
              union_indices=flex.union(len(observed), iselections=sample[crystal_model])
              union_observed = observed.select(union_indices)
              print ('done taking unions')
              # First index the union set with the central crystal model of the cluster
              self.known_crystal_models = None #[known_crystal_models[crystal_model]]
              from cctbx import crystal
              imagesets = datablock.extract_imagesets()
              explist_mean = ExperimentList()
              for i,imageset in enumerate(imagesets):
                exp = Experiment(imageset=imageset, 
                                 beam=imageset.get_beam(),
                                 detector=imageset.get_detector(), 
                                 goniometer=imageset.get_goniometer(),
                                 scan=imageset.get_scan(),
                                 crystal=known_crystal_models[crystal_model])
                explist_mean.append(exp)

              from exafel_project.ADSE13_25.indexing.indexer_iota import iota_indexer
              reidxr = iota_indexer(union_observed, imagesets,params=self.params)
              reidxr.calculate_fractional_hkl_from_Ainverse_q(reidxr.reflections, explist_mean)
              experiments_mean = explist_mean
              indexed_mean = reidxr.reflections
              indexed_mean['id'].set_selected(flex.size_t(range(len(indexed_mean))), crystal_model)
              print ('finished evaluating mean indexing results for crystal model ',crystal_model)
              # Now index with each each experimental model for each of the unioned observations
              dh_list = flex.double()
              failed_model_counter = 0
              for obs in all_experimental_models[crystal_model]:
                try:
                  explist = ExperimentList()
                  self.known_crystal_models = None #[obs.crystals()[0]]

                  # Make sure the crystal is rotated using the best_similarity_transformation
                  # with respect to the mean model. Otherwise dh values will be junk
                  from cctbx_orientation_ext import crystal_orientation
                  cryst_ref_ori = crystal_orientation(explist_mean.crystals()[0].get_A(), True)
                  cryst_tmp_ori = crystal_orientation(obs.crystals()[0].get_A(), True)
                  best_similarity_transform = cryst_tmp_ori.best_similarity_transformation(
                    other = cryst_ref_ori, fractional_length_tolerance = 10.00,
                    unimodular_generator_range=1)
                  cryst_tmp_ori_best=cryst_tmp_ori.change_basis(best_similarity_transform)
                  obs.crystals()[0].set_A(cryst_tmp_ori_best.reciprocal_matrix())

                  exp = Experiment(imageset=imageset,
                                   beam=imageset.get_beam(),
                                   detector=imageset.get_detector(),
                                   goniometer=imageset.get_goniometer(),
                                   scan=imageset.get_scan(),
                                   crystal=obs.crystals()[0])
                  explist.append(exp)
                  reidxr = iota_indexer(union_observed, imagesets,params=self.params)
                  reidxr.calculate_fractional_hkl_from_Ainverse_q(reidxr.reflections, explist)
                  experiments_tmp = explist
                  indexed_tmp = reidxr.reflections

                  # find dh = |h_frac - h_mean|
                  indexed_idxlist = [idx for idx,elem in enumerate(indexed_tmp['xyzobs.mm.value']) 
                                     if elem in indexed_mean['xyzobs.mm.value']]
                  dh_list_tmp = flex.double()
                  for idx in indexed_idxlist:
                    mean_list_idx = list(indexed_mean['xyzobs.mm.value']).index(indexed_tmp['xyzobs.mm.value'][idx])
                    x = indexed_mean['miller_index'][mean_list_idx]
                    y = indexed_tmp['fractional_miller_index'][idx]
                    dh = flex.double((x[0]-y[0],x[1]-y[1],x[2]-y[2])).norm()
                    if x == (0,0,0): continue

                    dh_list_tmp.append(dh)
                  dh_list.extend(dh_list_tmp)
                  print ('finished evaluating dh_list for crystal model ',crystal_model)
                except Exception as e:
                  print ('Reindexing with candidate lattices on union set failed')  
              # Get a sense of the variability in dh. Assign Z-score cutoff from there
              # try
              try:
                dh_mean_and_variance = flex.mean_and_variance(dh_list)
                dh_mean = dh_mean_and_variance.mean()
                dh_stdev = dh_mean_and_variance.unweighted_sample_standard_deviation()
                Z_cutoff = self.params.iota.random_sub_sampling.Z_cutoff
                # FIXME this part needs some though. What metric to use ?
                # Make sure cutoff is less than 0.5, else set it to 0.5
                #dh_cutoff = min(dh_mean+Z_cutoff*dh_stdev, 0.3)
                dh_cutoff = flex.median(dh_list)
                print ('dh_cutoff = ',dh_cutoff)
                # Now go through the spots indexed by the cluster center and reject if dh greater than Z_cutoff
                indexed_spots_idx = []
                for ii,refl in enumerate(indexed_mean):
                  dh = flex.double([refl['miller_index'][0]-refl['fractional_miller_index'][0],
                                refl['miller_index'][1]-refl['fractional_miller_index'][1],
                                refl['miller_index'][2]-refl['fractional_miller_index'][2]]).norm()
                  if dh < dh_cutoff and refl['miller_index'] != (0,0,0):
                    indexed_spots_idx.append(ii)
                indexed.extend(indexed_mean.select(flex.size_t(indexed_spots_idx)))
                # Need to append properly
                for iexpt,expt in enumerate(experiments_mean):
                  print ('APPENDING EXPERIMENT = ',crystal_model,iexpt)
                  experiments.append(expt)

              except:
                print ('dh_list calculation and outlier rejection failed')

            # Make sure crytal model numbers are in sequence, example 0,1,2 instead of 0,2,3 
            # when model 1 was not used for consensus part. Otherwise refine won't work
            max_id = flex.max(indexed['id'])
            original_ids = []
            for iid in range(0,max_id+1):
              if len(indexed.select(indexed['id']==iid)) != 0:
                original_ids.append(iid)

            for ii,iid in enumerate(original_ids):
              indexed['id'].set_selected(indexed['id'] == iid,ii)
            # Set back whatever PHIL parameter was supplied by user for outlier rejection and refinement
            self.params.indexing.stills.candidate_outlier_rejection=outlier_rejection_flag
            self.params.indexing.stills.refine_all_candidates=refine_all_candidates_flag

          try:
            experiments,indexed = self.refine(experiments, indexed)
          except Exception as e:
            print("Error refining", tag, str(e))
            if not self.params.dispatch.squash_errors: raise
            return
          try:
            if self.params.dispatch.integrate:
              integrated = self.integrate(experiments, indexed)
            else:
              print("Integration turned off. Exiting")
              return
          except Exception as e:
            print("Error integrating", tag, str(e))
            if self.params.dispatch.squash_errors: raise
            return
      else:
        print("IOTA based Indexing turned off. Exiting")
        return
    except Exception as e:
      print("Couldnt index using IOTA ", tag, str(e))
      if not self.params.dispatch.squash_errors: raise
      return

  def index(self, datablock, reflections):
    from dials.algorithms.indexing.indexer import indexer_base
    from time import time
    import copy
    st = time()

    logger.info('*' * 80)
    logger.info('Indexing Strong Spots')
    logger.info('*' * 80)

    imagesets = datablock.extract_imagesets()

    params = copy.deepcopy(self.params)
    # don't do scan-varying refinement during indexing
    params.refinement.parameterisation.scan_varying = False

    if hasattr(self, 'known_crystal_models'):
      known_crystal_models = self.known_crystal_models
    else:
      known_crystal_models = None

    if params.indexing.stills.method_list is None:
      idxr = indexer_base.from_parameters(
        reflections, imagesets, known_crystal_models=known_crystal_models,
        params=params)
      idxr.index()
    else:
      indexing_error = None
      for method in params.indexing.stills.method_list:
        params.indexing.method = method
        try:
          idxr = indexer_base.from_parameters(
            reflections, imagesets,
            params=params)
          idxr.index()
        except Exception as e:
          logger.info("Couldn't index using method %s"%method)
          if indexing_error is None:
            if e is None:
              e = Exception("Couldn't index using method %s"%method)
            indexing_error = e
        else:
          indexing_error = None
          break
      if indexing_error is not None:
        raise indexing_error

    indexed = idxr.refined_reflections
    experiments = idxr.refined_experiments

    if known_crystal_models is not None:
      from dials.array_family import flex
      filtered = flex.reflection_table()
      for idx in set(indexed['miller_index']):
        sel = indexed['miller_index'] == idx
        if sel.count(True) == 1:
          filtered.extend(indexed.select(sel))
      logger.info("Filtered duplicate reflections, %d out of %d remaining"%(len(filtered),len(indexed)))
      print("Filtered duplicate reflections, %d out of %d remaining"%(len(filtered),len(indexed)))
      indexed = filtered

    logger.info('')
    logger.info('Time Taken = %f seconds' % (time() - st))
    return experiments, indexed


  def index_with_iota(self, datablock, reflections):
    from exafel_project.ADSE13_25.indexing.indexer_iota import iota_indexer
    from time import time
    import copy
    st = time()

    logger.info('*' * 80)
    logger.info('Indexing Strong Spots')
    logger.info('*' * 80)

    imagesets = datablock.extract_imagesets()

    params = copy.deepcopy(self.params)
    # don't do scan-varying refinement during indexing
    params.refinement.parameterisation.scan_varying = False

    if hasattr(self, 'known_crystal_models'):
      known_crystal_models = self.known_crystal_models
    else:
      known_crystal_models = None

    if params.indexing.stills.method_list is None:
      idxr = iota_indexer.from_parameters(
        reflections, imagesets, known_crystal_models=known_crystal_models,
        params=params)
      idxr.index()
    else:
      indexing_error = None
      for method in params.indexing.stills.method_list:
        params.indexing.method = method
        try:
          idxr = iota_indexer.from_parameters(
            reflections, imagesets,
            params=params)
          idxr.index()
        except Exception as e:
          logger.info("Couldn't index using method %s"%method)
          if indexing_error is None:
            if e is None:
              e = Exception("Couldn't index using method %s"%method)
            indexing_error = e
        else:
          indexing_error = None
          break
      if indexing_error is not None:
        raise indexing_error

    indexed = idxr.reflections
    experiments = idxr.experiments

    if known_crystal_models is not None:
      from dials.array_family import flex
      filtered = flex.reflection_table()
      for idx in set(indexed['miller_index']):
        sel = indexed['miller_index'] == idx
        if sel.count(True) == 1:
          filtered.extend(indexed.select(sel))
      logger.info("Filtered duplicate reflections, %d out of %d remaining"%(len(filtered),len(indexed)))
      print("Filtered duplicate reflections, %d out of %d remaining"%(len(filtered),len(indexed)))
      indexed = filtered

    logger.info('')
    logger.info('Time Taken = %f seconds' % (time() - st))
    return experiments, indexed


if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script_iota()
    script.run()
  except Exception as e:
    halraiser(e)
