from __future__ import absolute_import, division
from __future__ import print_function
import math
import logging
logger = logging.getLogger(__name__)

from dials.util import log

debug_handle = log.debug_handle(logger)
info_handle = log.info_handle(logger)

import libtbx
from libtbx.utils import Sorry

from dials.algorithms.indexing.indexer import indexer_base
from dials.algorithms.indexing.known_orientation import indexer_known_orientation
from dials.algorithms.indexing.real_space_grid_search import indexer_real_space_grid_search
from dials.algorithms.indexing.fft3d import indexer_fft3d
from dials.algorithms.indexing.fft1d import indexer_fft1d

import iotbx.phil
from scitbx import matrix

from dials.array_family import flex
from cctbx import crystal, sgtbx, xray

from dxtbx.model import Crystal
from dxtbx.model.experiment_list import Experiment, ExperimentList

from dials.algorithms.indexing.indexer import max_cell_phil_str, index_only_phil_str,master_params 
from dials.algorithms.indexing.stills_indexer import stills_indexer

class iota_indexer(stills_indexer)
 
  def __init__(self, reflections, imagesets, params=None):
    '''Init function for iota_indexer is different from indexer_base in that
       _setup_symmetry function is not called. All features only work for stills'''
    stills_indexer.__init__(self, reflections, imagesets, params)
    self.reflections = reflections
    self.imagesets = imagesets
    if params is None: params = master_params
    self.params = params.indexing
    self.all_params = params
    self.refined_experiments = None
    self.hkl_offset = None
    
    if self.params.refinement_protocol.n_macro_cycles in ('auto', libtbx.Auto):
      self.params.refinement_protocol.n_macro_cycles = 1

    for imageset in imagesets[1:]:
      if imageset.get_detector().is_similar_to(self.imagesets[0].get_detector()):
        imageset.set_detector(self.imagesets[0].get_detector())

    if 'flags' in self.reflections:
      strong_sel = self.reflections.get_flags(self.reflections.flags.strong)
      if strong_sel.count(True) > 0:
        self.reflections = self.reflections.select(strong_sel)
    if 'flags' not in self.reflections or strong_sel.count(True) == 0:
      # backwards compatibility for testing
      self.reflections.set_flags(
        flex.size_t_range(len(self.reflections)), self.reflections.flags.strong)

    #self._setup_symmetry()
    self.d_min = None
    self.setup_indexing()

  @staticmethod
  def from_parameters(reflections, imagesets,
                      known_crystal_models=None, params=None):
    '''Sets up indexer object that will be used for indexing '''
    if params is None:
      params = master_params
    
    if params.indexing.basis_vector_combinations.max_refine is libtbx.Auto:
      params.indexing.basis_vector_combinations.max_refine = 5

    # Ensure the indexer and downstream applications treat this as set of stills
    from dxtbx.imageset import ImageSet 
    reset_sets = []
    for i in range(len(imagesets)):
      imagesweep = imagesets.pop(0)
      imageset = ImageSet(imagesweep.data(), imagesweep.indices())
      imageset.set_scan(None)
      imageset.set_goniometer(None)
      reset_sets.append(imagesets)
    imagesets.extend(reset_sets)

    if known_crystal_models is not None:
      from dials.algorithms.indexing.known_orientation \
        import indexer_known_orientation
      idxr = indexer_known_orientation(
        reflections, imagesets, params, known_crystal_models)
    elif params.indexing.method == "fft3d":
      idxr = iota_indexer_fft3d(reflections, imagesets, params=params)
    elif params.indexing.method == "fft1d":
      idxr = iota_indexer_fft1d(reflections, imagesets, params=params)
    elif params.indexing.method == "real_space_grid_search":
      idxr = iota_indexer_real_space_grid_search(reflections, imagesets, params=params)
    return idxr

  def index(self):
    ''' This step calls find_candidate_orientation_matrices, best_orientation_matrix,
        apply_symmetry. Eventually returns all indexing solutions'''
    if self.params.refinement.protocol.n_macro_cycles > 1:
      raise Sorry('For Stills, set refinement_protocol.n_macro_cycles = 1')

    experiments = ExperimentList()
    had_refinement_error = False
    have_similar_crystal_models = False

    while True:
      self.d_min = self.params.refinement_protocol.d_min_start
      if had_refinement_error or have_similar_crystal_models:
        break

      if max_lattices is not None and len(experiments) >= max_lattices:
        break

      if len(experiments) > 0:
        cutoff_fraction

  def index_reflections(self, experiments, reflections):
    ''' Does actual indexing'''
    pass

class iota_indexer_fft1d(iota_indexer, indexer_fft1d):
  '''Mixin class with fft1d and iota_indexer way of identifying indexing solutions'''
  pass

class iota_indexer_real_space_grid_search(iota_indexer, indexer_real_space_grid_search):
  '''Mixin class with real_space_grid_search and iota_indexer way of identifying indexing solutions'''
  pass

class iota_indexer_fft3d(iota_indexer, indexer_fft3d):
  '''Mixin class with fft3d and iota_indexer way of identifying indexing solutions'''
  pass
