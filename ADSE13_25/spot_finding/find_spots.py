from __future__ import absolute_import, division, print_function

import logging
import os
logger = logging.getLogger('exafel.find_spots')
from dxtbx.datablock import DataBlockFactory
from libtbx.utils import Abort, Sorry
from libtbx.phil import parse


# Imports for LS49
from dials.command_line.stills_process import Script, Processor, control_phil_str, dials_phil_str, program_defaults_phil_str

phil_scope = parse(control_phil_str + dials_phil_str, process_includes=True).fetch(parse(program_defaults_phil_str))


message = ''' This is a specialized version of stills_process which only does spot_finding. Difference is it only
              prints out the datablock/strong pickle files of those images which have hits on them
'''


def do_import(filename):
  logger.info("Loading %s"%os.path.basename(filename))
  datablocks = DataBlockFactory.from_filenames([filename])
  if len(datablocks) == 0:
    try:
      datablocks = DataBlockFactory.from_json_file(filename)
    except ValueError:
      raise Abort("Could not load %s"%filename)

  if len(datablocks) == 0:
    raise Abort("Could not load %s"%filename)
  if len(datablocks) > 1:
    raise Abort("Got multiple datablocks from file %s"%filename)

  # Ensure the indexer and downstream applications treat this as set of stills
  reset_sets = []

  from dxtbx.imageset import ImageSetFactory
  for imageset in datablocks[0].extract_imagesets():
    imageset = ImageSetFactory.imageset_from_anyset(imageset)
    imageset.set_scan(None)
    imageset.set_goniometer(None)
    reset_sets.append(imageset)

  return DataBlockFactory.from_imageset(reset_sets)[0]


class SpotFinding_Script(Script):
  def run(self):
    '''Execute the script.'''
    from dials.util import log
    from time import time
    import copy

    # Parse the command line
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
      info='exafel.spot_finding.log',
      debug='exafel.spot_finding.debug.log')

    bad_phils = [f for f in all_paths if os.path.splitext(f)[1] == ".phil"]
    if len(bad_phils) > 0:
      self.parser.print_help()
      logger.error('Error: the following phil files were not understood: %s'%(", ".join(bad_phils)))
      return

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
    if True: #pre_import:
      # Handle still imagesets by breaking them apart into multiple datablocks
      # Further handle single file still imagesets (like HDF5) by tagging each
      # frame using its index

      datablocks = [do_import(path) for path in all_paths]


      indices = []
      basenames = []
      datablock_references = []
      for datablock in datablocks:
        for j, imageset in enumerate(datablock.extract_imagesets()):
          paths = imageset.paths()
          for i in xrange(len(imageset)):
            datablock_references.append(datablock)
            indices.append((j,i))
            basenames.append(os.path.splitext(os.path.basename(paths[i]))[0])
      tags = []
      frame_number = []
      for (j,i), basename in zip(indices, basenames):
        if basenames.count(basename) > 1:
          if len(set([idx[0] for idx in indices if idx[0]==j])) > 1:
            tags.append("%s_%05d_%05d"%(basename, j, i))
          else:
            tags.append("%s_%05d"%(basename, i))
            frame_number.append(i)
        else:
          tags.append(basename)

      # Wrapper function
      def do_work(i, item_list):
        processor = SpotFinding_Processor(copy.deepcopy(params), composite_tag = "%04d"%i, rank = i)
        from dials.array_family import flex
        all_spots_from_rank = flex.reflection_table()

        for item in item_list:
          tag, (imgset_id, img_id), datablock = item
          #from IPython import embed; embed(); exit()
          imageset = datablock.extract_imagesets()[imgset_id]
          subset = imageset[img_id:img_id+1]
          try:
            update_geometry(subset)
          except RuntimeError as e:
            logger.warning("Error updating geometry on item %s, %s"%(str(tag), str(e)))
            continue

          if self.reference_detector is not None:
            from dxtbx.model import Detector
            subset.set_detector(
              Detector.from_dict(self.reference_detector.to_dict()),
              index=0)
          datablock = DataBlockFactory.from_imageset(subset)[0]

          try:
            refl_table = processor.process_datablock(tag, datablock, img_id)
            if refl_table is not None:
              all_spots_from_rank.extend(refl_table)
            #all_spots_from_rank.extend(processor.process_datablock(tag, datablock, img_id))
          except Exception as e:
            logger.warning("Unhandled error on item %s, %s"%(str(tag), str(e)))
        processor.finalize()
        #from IPython import embed; embed(); exit()
        return all_spots_from_rank

      iterable = zip(tags, indices, datablock_references)

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
      log.config(params.verbosity, info=info_path, debug=debug_path)

      subset = [item for i, item in enumerate(iterable) if (i+rank)%size == 0]
      all_spots_from_rank = do_work(rank, subset)
      all_spots_rank0 = comm.gather(all_spots_from_rank, root=0)
      from dials.array_family import flex
      all_spots = flex.reflection_table()
      for ii,refl_table in enumerate(all_spots_rank0):
        all_spots.extend(refl_table)

      if rank == 0:
        from dials.algorithms.spot_finding import per_image_analysis
        from cStringIO import StringIO
        s = StringIO()
        # Assuming one datablock. Might be dangerous
        # FIXME
        for i, imageset in enumerate(datablocks[0].extract_imagesets()):
          print("Number of centroids per image for imageset %i:" %i, file=s)
          #from IPython import embed; embed(); exit()
          stats = custom_stats_imageset(
            imageset, all_spots.select(all_spots['id'] == i))
          per_image_analysis.print_table(stats)
        logger.info(s.getvalue())



class SpotFinding_Processor(Processor):
  def process_datablock(self, tag, datablock,img_id):
    import os
    if not self.params.output.composite_output:
      self.setup_filenames(tag)
    self.tag = tag
    self.img_id = img_id
    self.debug_start(tag)

#    if not self.params.output.composite_output and self.params.output.datablock_filename:
#      from dxtbx.datablock import DataBlockDumper
#      dump = DataBlockDumper(datablock)
#      dump.as_json(self.params.output.datablock_filename)

    # Do spotfinding
    try:
      self.debug_write("spotfind_start")
      observed = self.find_spots(datablock)
      return observed
    except Exception as e:
      print("Error spotfinding", tag, str(e))
      return None

  def find_spots(self, datablock):
    from time import time
    from dials.array_family import flex
    st = time()

    logger.info('*' * 80)
    logger.info('Finding Strong Spots')
    logger.info('*' * 80)

    # Find the strong spots
    observed = flex.reflection_table.from_observations(datablock, self.params)

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
        from dxtbx.datablock import DataBlockDumper
        dump = DataBlockDumper(datablock)
        dump.as_json(self.params.output.datablock_filename)

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
      reflections.select(reflections['img_id']==i+start), i=i+start,
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
  from dials.util import halraiser
  try:
    script = SpotFinding_Script()
    script.run()
  except Exception as e:
    halraiser(e)
