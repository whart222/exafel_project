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
          # No outlier rejection or refinement should be done for the candidate basis vectors
          outlier_rejection_flag=self.params.indexing.stills.candidate_outlier_rejection
          refine_all_candidates_flag=self.params.indexing.stills.refine_all_candidates
          if self.params.iota.random_sub_sampling.no_outlier_rejection_and_candidates_refinement:
            self.params.indexing.stills.candidate_outlier_rejection=False
            self.params.indexing.stills.refine_all_candidates=False
          for trial in range(self.params.iota.random_sub_sampling.ntrials):
            flex.set_random_seed(trial+1001)
            observed_sample = observed.select(flex.random_selection(len(observed), int(len(observed)*self.params.iota.random_sub_sampling.fraction_sub_sample)))
            try:
              print ('IOTA: SUM_INTENSITY_VALUE',sum(observed_sample['intensity.sum.value']), ' ',trial)
              experiments_tmp, indexed_tmp = self.index(datablock, observed_sample)
              experiments_list.append(experiments_tmp)
            except:
              print('Indexing failed for some reason')
          if self.params.iota.random_sub_sampling.consensus_function == 'unit_cell':
            from exafel_project.ADSE13_25.consensus_functions import get_uc_consensus as get_consensus
            known_crystal_models = get_consensus(experiments_list, show_plot=self.params.iota.random_sub_sampling.show_plot, return_only_first_indexed_model=True)
          self.known_crystal_models = known_crystal_models
          print ('IOTA: Reindexing with best chosen crystal model')
          # Set back whatever PHIL parameter was supplied by user for outlier rejection and refinement
          self.params.indexing.stills.candidate_outlier_rejection=outlier_rejection_flag
          self.params.indexing.stills.refine_all_candidates=refine_all_candidates_flag
          experiments, indexed = self.index(datablock, observed)
          print('fraction subsampled = %5.2f with %d indexed spots ' %(self.params.iota.random_sub_sampling.fraction_sub_sample,len(indexed)))
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


if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script_iota()
    script.run()
  except Exception as e:
    halraiser(e)
