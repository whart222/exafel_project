from __future__ import absolute_import, division
from __future__ import print_function
import logging
logger = logging.getLogger(__name__)

from dials.util import log

import libtbx
from libtbx.utils import Sorry

#from dials.algorithms.indexing.indexer import indexer_base
#from dials.algorithms.indexing.known_orientation import indexer_known_orientation
#from dials.algorithms.indexing.real_space_grid_search import indexer_real_space_grid_search
#from exafel_project.ADSE13_25.indexing.real_space_grid_smart_search import indexer_real_space_grid_smart_search
#from dials.algorithms.indexing.fft3d import indexer_fft3d
#from dials.algorithms.indexing.fft1d import indexer_fft1d
from dials_algorithms_indexing_ext import *
from dials.algorithms.indexing.symmetry import SymmetryHandler

#Aug_Refactor
from dials.algorithms.indexing.known_orientation import  IndexerKnownOrientation
from dials.algorithms.indexing.lattice_search import BasisVectorSearch, LatticeSearch

import iotbx.phil # implicit import
from scitbx import matrix

from dials.array_family import flex
from cctbx import crystal

from dxtbx.model import Crystal
from dxtbx.model.experiment_list import ExperimentList

#from dials.algorithms.indexing.indexer import master_params
from dials.algorithms.indexing.stills_indexer import StillsIndexer
import pkg_resources

class iota_indexer(StillsIndexer):

  def __init__(self, reflections, experiments, params=None):
    '''Init function for iota_indexer is different from indexer_base in that
       _setup_symmetry function is not called. All features only work for stills'''


    # FIXME this should not be called the stills_indexer __init__ method
    # FIXME need to write own __init__ function
    #stills_indexer.__init__(self, reflections, imagesets, params)
    self.reflections = reflections
    self.experiments = experiments
    #if params is None: params = master_params
    self.params = params.indexing
    self.all_params = params
    self.refined_experiments = None
    self.hkl_offset = None
    #Aug_Refactor
    #if self.params.index_assignment.method == "local":
    #  self._assign_indices = assign_indices.AssignIndicesLocal(
    #    epsilon=self.params.index_assignment.local.epsilon,
    #    delta=self.params.index_assignment.local.delta,
    #    l_min=self.params.index_assignment.local.l_min,
    #    nearest_neighbours=self.params.index_assignment.local.nearest_neighbours)
    #else:
    #  self._assign_indices = assign_indices.AssignIndicesGlobal(
    #    tolerance=self.params.index_assignment.simple.hkl_tolerance)
    ###

    if self.params.refinement_protocol.n_macro_cycles in ('auto', libtbx.Auto):
      self.params.refinement_protocol.n_macro_cycles = 1

    for expt in self.experiments[1:]:
      if expt.detector.is_similar_to(self.experiments[0].detector):
        expt.detector=self.experiments[0].detector

    if 'flags' in self.reflections:
      strong_sel = self.reflections.get_flags(self.reflections.flags.strong)
      if strong_sel.count(True) > 0:
        self.reflections = self.reflections.select(strong_sel)
    if 'flags' not in self.reflections or strong_sel.count(True) == 0:
      # backwards compatibility for testing
      self.reflections.set_flags(
        flex.size_t_range(len(self.reflections)), self.reflections.flags.strong)

    self._setup_symmetry()
    self.d_min = None
    self.setup_indexing()

  @staticmethod
  def from_parameters(reflections, experiments,
                      known_crystal_models=None, params=None):
    '''Sets up indexer object that will be used for indexing '''
#    if params is None:
#      params = master_params

    if known_crystal_models is not None:
      #from dials.algorithms.indexing.known_orientation \
      #     import indexer_known_orientation
      #idxr = iota_indexer_known_orientation(
      #  reflections, experiments, params, known_crystal_models)

      idxr=IOTA_StillsIndexerKnownOrientation(reflections, experiments, params, known_crystal_models)

      #idxr = indexer_known_orientation(
      #  reflections, imagesets, params, known_crystal_models)
    else:
      has_stills = False
      has_sweeps = False
      for expt in experiments:
        if expt.goniometer is None or expt.scan is None or \
            expt.scan.get_oscillation()[1] == 0:
          if has_sweeps:
            raise Sorry("Please provide only stills or only sweeps, not both")
          has_stills = True
        else:
          if has_stills:
            raise Sorry("Please provide only stills or only sweeps, not both")
          has_sweeps = True
      assert not (has_stills and has_sweeps)
      use_stills_indexer = has_stills

      if not (params.indexing.stills.indexer is libtbx.Auto or params.indexing.stills.indexer.lower() == 'auto'):
        if params.indexing.stills.indexer == 'stills':
          use_stills_indexer = True
        elif params.indexing.stills.indexer == 'sweeps':
          use_stills_indexer = False
        else:
          assert False


    if params.indexing.basis_vector_combinations.max_refine is libtbx.Auto:
      params.indexing.basis_vector_combinations.max_refine = 5

    # Ensure the indexer and downstream applications treat this as set of stills
    from dxtbx.imageset import ImageSet
    # DIALS 2.0 stuff here
    for experiment in experiments:
      experiment.imageset=ImageSet(experiment.imageset.data(), experiment.imageset.indices())
      experiment.imageset.set_scan(None)
      experiment.imageset.set_goniometer(None)
      experiment.scan=None
      experiment.goniometer=None

    # Old code prior to DIALS 2.0
    #reset_sets = []
    #for i in range(len(imagesets)):
    #  imagesweep = imagesets.pop(0)
    #  imageset = ImageSet(imagesweep.data(), imagesweep.indices())
    #  imageset.set_scan(None)
    #  imageset.set_goniometer(None)
    #  reset_sets.append(imageset)
    #imagesets.extend(reset_sets)

    if known_crystal_models is not None:
      print ('Known Orientation Indexing')
      return idxr
      #from dials.algorithms.indexing.known_orientation \
      #  import indexer_known_orientation
      #idxr = indexer_known_orientation(
      #  reflections, experiments, params, known_crystal_models)
    for entry_point in pkg_resources.iter_entry_points("dials.index.basis_vector_search"):
      if params.indexing.method==entry_point.name:
        idxr=IOTA_StillsIndexerBasisVectorSearch(reflections, experiments, params=params)
        return idxr
          
    #elif params.indexing.method == "fft3d":
    #  idxr = iota_indexer_fft3d(reflections, experiments, params=params)
    #elif params.indexing.method == "fft1d":
    #  idxr = iota_indexer_fft1d(reflections, experiments, params=params)
    #elif params.indexing.method == "real_space_grid_search":
    #  idxr = iota_indexer_real_space_grid_search(reflections, experiments, params=params)
    #elif params.indexing.method == "real_space_grid_smart_search":
    #  idxr = iota_indexer_real_space_grid_smart_search(reflections, experiments, params=params)
    #return idxr

  def index(self, provided_experiments=None,debug=False):
    ''' This step does  1. find_lattices (via a method like fft1d)
                        2. Assign hkl indices (through index_reflections)
                        3. Housekeeping like apply_symmetry, discard too similar models
    '''

    experiments = ExperimentList()
    have_similar_crystal_models = False

    self.d_min = self.params.refinement_protocol.d_min_start

    # Find lattices i.e the basis vectors & unit cell params
    # Possible to index multiple lattices  ??
    while True:
      max_lattices = self.params.multiple_lattice_search.max_lattices
      if max_lattices is not None and len(experiments) >= max_lattices:
        break 
      n_lattices_previous_cycle = len(experiments)
      if len(experiments) == 0:
        experiments.extend(self.find_lattices())
      else:
        try:
          new = self.find_lattices()
          experiments.extend(new)
        except Sorry:
          print ('Indexing remaining reflections failed')
          break


      if len(experiments) == 0:
        raise Sorry("No suitable lattice could be found.")

      # Initialize id values as -1 since no indexing has been done yet
      self.reflections['id'] = flex.int(len(self.reflections), -1)

      # Now index reflections
      self.index_reflections(experiments, self.reflections, debug=debug)

      # Housekeeping. Apply symmetry
      self._apply_symmetry_post_indexing(experiments, self.reflections, n_lattices_previous_cycle) 

      # Aug_Refactor :: probably unnecessary to remove stuff below but still doing so to adapt to new style
      # Never mind, might keep it to have IOTA stuff working
      '''
      target_space_group = self.target_symmetry_primitive.space_group()
      for i_cryst, cryst in enumerate(experiments.crystals()):
        if i_cryst >= n_lattices_previous_cycle:
          new_cryst, cb_op_to_primitive = self.apply_symmetry(
                                        cryst, target_space_group)
          if provided_experiments is None:
            if self.cb_op_primitive_inp is not None:
              new_cryst = new_cryst.change_basis(self.cb_op_primitive_inp)
              logger.info(new_cryst.get_space_group().info())
            cryst.update(new_cryst)
            cryst.set_space_group(
                self.params.known_symmetry.space_group.group())
          for i_expt, expt in enumerate(experiments):
            if expt.crystal is not cryst:
              continue
            if not cb_op_to_primitive.is_identity_op():
              miller_indices = self.reflections['miller_index'].select(
                  self.reflections['id'] == i_expt)

              if provided_experiments is None:
                miller_indices = cb_op_to_primitive.apply(miller_indices)
              self.reflections['miller_index'].set_selected(
                  self.reflections['id'] == i_expt, miller_indices)

            if self.cb_op_primitive_inp is not None:
              miller_indices = self.reflections['miller_index'].select(
                  self.reflections['id'] == i_expt)

              if provided_experiments is None:
                miller_indices = self.cb_op_primitive_inp.apply(miller_indices)
              self.reflections['miller_index'].set_selected(
                  self.reflections['id'] == i_expt, miller_indices)
              # IOTA
              from scitbx.matrix import sqr
              hklfrac=flex.mat3_double(len(miller_indices), sqr(cryst.get_A()).inverse())*self.reflections['rlp'].select(self.reflections['id']==i_expt)
              self.reflections['fractional_miller_index'].set_selected(self.reflections['id']==i_expt, hklfrac)
    '''

    logger.info("\nIndexed crystal models:")  
    self.show_experiments(experiments, self.reflections, d_min=self.d_min)

    # Discard nearly overlapping lattices
    # difference_rotation_matrix_axis_angle function is there still in DIALS 2.0 so no need to change anything below 
    if len(experiments) > 1:
      from dials.algorithms.indexing.compare_orientation_matrices \
        import difference_rotation_matrix_axis_angle
      cryst_b = experiments.crystals()[-1]
      have_similar_crystal_models = False
      for i_a, cryst_a in enumerate(experiments.crystals()[:-1]):
        R_ab, axis, angle, cb_op_ab = \
        difference_rotation_matrix_axis_angle(cryst_a, cryst_b)
        min_angle = self.params.multiple_lattice_search.minimum_angular_separation
        if abs(angle) < min_angle: # degrees
          logger.info("Crystal models too similar, rejecting crystal %i:" %(
              len(experiments)))
          logger.info("Rotation matrix to transform crystal %i to crystal %i" %(
              i_a+1, len(experiments)))
          logger.info(R_ab)
          logger.info("Rotation of %.3f degrees" %angle + " about axis (%.3f, %.3f, %.3f)" %axis)
          have_similar_crystal_models = True
          del experiments[-1]
          break

    self.indexed_reflections = (self.reflections['id'] > -1)
    self.experiments = experiments


  def index_reflections(self, experiments, reflections,debug=False):
    ''' Assigns hkl values to reflections'''

    # Not needed below as I call assign_hkl_to_reflections
    #self._assign_indices(reflections, experiments, d_min=self.d_min)
    params_simple = self.params.index_assignment.simple
    self.assign_hkl_to_reflections(reflections, experiments, self.d_min,
                      tolerance = params_simple.hkl_tolerance, debug=debug)

    if self.hkl_offset is not None and self.hkl_offset != (0,0,0):
      reflections['miller_index'] = apply_hkl_offset(
        reflections['miller_index'], self.hkl_offset)
      self.hkl_offset = None

  def assign_hkl_to_reflections(self, reflections, experiments, d_min=None, tolerance=0.3,debug=False):
    ''' Function to assign hkl values to reflections on a shot. Uses underlying c++
      function AssignIndices() available in DIALS'''
    from cctbx.array_family import flex

    reciprocal_lattice_points = reflections['rlp']
    reflections['miller_index'] = flex.miller_index(len(reflections), (0,0,0))
    # IOTA
    reflections['fractional_miller_index'] = flex.vec3_double(len(reflections),(0.0,0.0,0.0))

    if d_min is not None:
      d_spacings = 1/reciprocal_lattice_points.norms()
      inside_resolution_limit = d_spacings > d_min
    else:
      inside_resolution_limit = flex.bool(reciprocal_lattice_points.size(), True)

    sel = inside_resolution_limit & (reflections['id'] == -1)
    isel = sel.iselection()
    rlps = reciprocal_lattice_points.select(isel)
    refs = reflections.select(isel)
    phi = refs['xyzobs.mm.value'].parts()[2]

    diffs = []
    norms = []
    hkl_ints = []

    UB_matrices = flex.mat3_double([cm.get_A() for cm in experiments.crystals()])
    imgset_ids = reflections['imageset_id'].select(sel)

    for i_imgset, imgset in enumerate(experiments.imagesets()):
      sel_imgset = (imgset_ids == i_imgset)
      result = AssignIndices(
        rlps.select(sel_imgset), phi.select(sel_imgset), UB_matrices, tolerance=tolerance)

      miller_indices = result.miller_indices()
      crystal_ids = result.crystal_ids()
      expt_ids = flex.int(crystal_ids.size(), -1)
      for i_cryst, cryst in enumerate(experiments.crystals()):
        sel_cryst = (crystal_ids == i_cryst)
        for i_expt in experiments.where(
          crystal=cryst, imageset=imgset):
          expt_ids.set_selected(sel_cryst, i_expt)

      reflections['miller_index'].set_selected(isel.select(sel_imgset), miller_indices)
      reflections['id'].set_selected(isel.select(sel_imgset), expt_ids)
      reflections.set_flags(
        reflections['miller_index'] != (0,0,0), reflections.flags.indexed)
      reflections['id'].set_selected(reflections['miller_index'] == (0,0,0), -1)

  def calculate_fractional_hkl_from_Ainverse_q(self, reflections, experiments,debug=False):
    ''' Calculate hkl_frac = A^-1*q. Will also calculate the integer hkl values '''
    assert len(experiments.crystals()) == 1, 'Should have only one crystal model'
    #self.map_centroids_to_reciprocal_space(observed, )
    from scitbx.matrix import sqr
    Ainverse = flex.mat3_double(len(reflections), sqr(experiments.crystals()[0].get_A()).inverse())
    q = reflections['rlp']
    hkl_frac = Ainverse*q
    hkl = hkl_frac.iround()
    reflections['miller_index'] = flex.miller_index(len(reflections),(0,0,0))
    reflections['fractional_miller_index'] = flex.vec3_double(len(reflections),(0.0,0.0,0.0))
    reflections['miller_index'] = flex.miller_index(list(hkl))
    reflections['fractional_miller_index'] = hkl_frac
    reflections.set_flags(reflections['miller_index'] != (0,0,0), reflections.flags.indexed)
    reflections['id'].set_selected(flex.size_t(range(len(reflections))), 0)
    # Add predicted reflections
    from dials.algorithms.refinement.prediction.managed_predictors import ExperimentsPredictorFactory
    ref_predictor = ExperimentsPredictorFactory.from_experiments(experiments, force_stills=experiments.all_stills())
    reflections = ref_predictor(reflections)
    reflections['id'].set_selected(flex.size_t(range(len(reflections))), -1)

    def _apply_symmetry_post_indexing(self, experiments, reflections, n_lattices_previous_cycle):
      # Taken from DIALS post Aug_Refactor
      # now apply the space group symmetry only after the first indexing
      # need to make sure that the symmetrized orientation is similar to the P1 model
      for cryst in experiments.crystals()[n_lattices_previous_cycle:]:
        new_cryst, cb_op_to_primitive = self._symmetry_handler.apply_symmetry(cryst)
        if self._symmetry_handler.cb_op_primitive_inp is not None:
            new_cryst = new_cryst.change_basis(self._symmetry_handler.cb_op_primitive_inp)
        cryst.update(new_cryst)
        cryst.set_space_group(self.params.known_symmetry.space_group.group())
        for i_expt, expt in enumerate(experiments):
            if expt.crystal is not cryst:
                continue
            if not cb_op_to_primitive.is_identity_op():
                miller_indices = reflections["miller_index"].select(reflections["id"] == i_expt)
                miller_indices = cb_op_to_primitive.apply(miller_indices)
                reflections["miller_index"].set_selected(reflections["id"] == i_expt, miller_indices)
            if self._symmetry_handler.cb_op_primitive_inp is not None:
                miller_indices = reflections["miller_index"].select(reflections["id"] == i_expt)
                miller_indices = self._symmetry_handler.cb_op_primitive_inp.apply(miller_indices)
                reflections["miller_index"].set_selected(reflections["id"] == i_expt, miller_indices)
                # IOTA
                from scitbx.matrix import sqr
                hklfrac=flex.mat3_double(len(miller_indices), sqr(cryst.get_A()).inverse())*self.reflections['rlp'].select(self.reflections['id']==i_expt)
                self.reflections['fractional_miller_index'].set_selected(self.reflections['id']==i_expt, hklfrac)


class IOTA_StillsIndexerKnownOrientation(IndexerKnownOrientation, iota_indexer):
    pass


class IOTA_StillsIndexerBasisVectorSearch(iota_indexer, BasisVectorSearch):
    def __init__(self, reflections, experiments, params):
      BasisVectorSearch.__init__(self,reflections, experiments, params)
      iota_indexer.__init__(self,reflections, experiments, params)
    #pass


class IOTA_StillsIndexerLatticeSearch(iota_indexer, LatticeSearch):
    pass

## Post refactor
#class iota_indexer_fft1d(iota_indexer, indexer_fft1d):
#  '''Mixin class with fft1d and iota_indexer way of identifying indexing solutions'''
#  pass
#
#class iota_indexer_real_space_grid_search(iota_indexer, indexer_real_space_grid_search):
#  '''Mixin class with real_space_grid_search and iota_indexer way of identifying indexing solutions'''
#  pass
#
#class iota_indexer_real_space_grid_smart_search(iota_indexer, indexer_real_space_grid_smart_search):
#  '''Mixin class with real_space_grid_smart_search and iota_indexer way of identifying indexing solutions'''
#  pass
#
#class iota_indexer_fft3d(iota_indexer, indexer_fft3d):
#  '''Mixin class with fft3d and iota_indexer way of identifying indexing solutions'''
#  pass
#
#class iota_indexer_known_orientation(indexer_known_orientation, iota_indexer):
#  pass
