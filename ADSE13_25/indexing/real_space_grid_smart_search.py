# Trying to be a bit smarter about the real_space_grid_search in DIALS by first doing a coarse grid search
# with the default of 0.02 grid size and then gradually reducing the grid size to 0.002 or smaller ?
# Looking to reduce the number of calls to compute functionals which is what eats up so much time

from __future__ import absolute_import, division, print_function
import math
import logging

logger = logging.getLogger(__name__)

from scitbx import matrix
from scitbx.array_family import flex
from dials.algorithms.indexing.indexer import indexer_base, optimise_basis_vectors
from dials.algorithms.indexing.indexer import is_approximate_integer_multiple
from dxtbx.model.experiment_list import Experiment, ExperimentList
from dials.algorithms.indexing.real_space_grid_search import indexer_real_space_grid_search

class indexer_real_space_grid_smart_search(indexer_real_space_grid_search):
  def __init__ (self, reflections, experiments, params):
    super(indexer_real_space_grid_smart_search, self).__init__(reflections, experiments, params)

  def find_lattices(self):
    #self.real_space_grid_search()
    self.real_space_grid_smart_search()
    crystal_models = self.candidate_crystal_models
    experiments = ExperimentList()
    for cm in crystal_models:
      for expt in self.experiments:
        experiments.append(Experiment(imageset=expt.imageset,
                        beam=expt.beam,
                        detector=expt.detector,
                        goniometer=expt.goniometer,
                        scan=expt.scan,
                        crystal=cm))
    return experiments

  def get_finegrained_SST(self, coarse_sampling_grid=0.02):
    d_min = self.params.refinement_protocol.d_min_start
    sel = self.reflections["id"] == -1
    if d_min is not None:
      sel &= 1 / self.reflections["rlp"].norms() > d_min
    reciprocal_lattice_points = self.reflections["rlp"].select(sel)
    print("Indexing from %i reflections COARSE" % len(reciprocal_lattice_points))
    def compute_functional(vector):
      two_pi_S_dot_v = 2 * math.pi * reciprocal_lattice_points.dot(vector)
      return flex.sum(flex.cos(two_pi_S_dot_v))

    from rstbx.array_family import flex
    from rstbx.dps_core import SimpleSamplerTool

    assert self.target_symmetry_primitive is not None
    assert self.target_symmetry_primitive.unit_cell() is not None
    SST = SimpleSamplerTool(coarse_sampling_grid)
    SST.construct_hemisphere_grid(SST.incr)
    cell_dimensions = self.target_symmetry_primitive.unit_cell().parameters()[:3]
    unique_cell_dimensions = set(cell_dimensions)
    print("Number of search vectors COARSE: %i"
                 %(len(SST.angles) * len(unique_cell_dimensions)))
    vectors = flex.vec3_double()
    function_values = flex.double()
    import time
    time1=time.time()
    SST_all_angles = flex.Direction()
    for i, direction in enumerate(SST.angles):
      for l in unique_cell_dimensions:
        v = matrix.col(direction.dvec) * l
        f = compute_functional(v.elems)
        vectors.append(v.elems)
        function_values.append(f)
        SST_all_angles.append(direction)
    time2=time.time()
    print ('COARSE GRID SEARCH TIME=',time2-time1)

    perm = flex.sort_permutation(function_values, reverse=True)
    vectors = vectors.select(perm)
    function_values = function_values.select(perm)

    unique_vectors = []
    unique_indices = []
    i = 0
    while len(unique_vectors) < 30:
      v = matrix.col(vectors[i])
      is_unique = True
      if i > 0:
        for v_u in unique_vectors:
          if v.length() < v_u.length():
            if is_approximate_integer_multiple(v, v_u):
              is_unique = False
              break
          elif is_approximate_integer_multiple(v_u, v):
            is_unique = False
            break
      if is_unique:
        unique_vectors.append(v)
        unique_indices.append(perm[i])
      i += 1

    # Evaluate which SST angles contributed to the unique vectors
    SST_filter = flex.Direction()
    for v in unique_indices:
      direction=SST.angles[v//len(unique_cell_dimensions)]
      SST_filter.append(direction)
    SST.construct_hemisphere_grid_finegrained(0.0005, coarse_sampling_grid, SST_filter)
    #from IPython import embed; embed(); exit()
    return SST

  def real_space_grid_smart_search(self):
    """ overloading the real_space_grid_search method to use a sparser SST with highly fine grid after initially doing
        a sparse grid to evaluate grid points of interest"""
    d_min = self.params.refinement_protocol.d_min_start
    sel = self.reflections["id"] == -1
    if d_min is not None:
      sel &= 1 / self.reflections["rlp"].norms() > d_min
    reciprocal_lattice_points = self.reflections["rlp"].select(sel)
    print("Indexing from %i reflections FINE " % len(reciprocal_lattice_points))
    def compute_functional(vector):
      two_pi_S_dot_v = 2 * math.pi * reciprocal_lattice_points.dot(vector)
      return flex.sum(flex.cos(two_pi_S_dot_v))

    from rstbx.array_family import flex
    from rstbx.dps_core import SimpleSamplerTool

    assert self.target_symmetry_primitive is not None
    assert self.target_symmetry_primitive.unit_cell() is not None

    SST = self.get_finegrained_SST()
    cell_dimensions = self.target_symmetry_primitive.unit_cell().parameters()[:3]
    unique_cell_dimensions = set(cell_dimensions)
    print("Number of search vectors FINE : %i"
                 %(len(SST.finegrained_angles) * len(unique_cell_dimensions)))
    vectors = flex.vec3_double()
    function_values = flex.double()
    import time
    time1=time.time()
    for i, direction in enumerate(SST.finegrained_angles):
      for l in unique_cell_dimensions:
        v = matrix.col(direction.dvec) * l
        f = compute_functional(v.elems)
        vectors.append(v.elems)
        function_values.append(f)
    time2=time.time()
    print ('FINE GRID SEARCH TIME=',time2-time1)
    #import pdb; pdb.set_trace()

    perm = flex.sort_permutation(function_values, reverse=True)
    vectors = vectors.select(perm)
    function_values = function_values.select(perm)

    unique_vectors = []
    unique_indices = []
    i = 0
    while len(unique_vectors) < 30:
      v = matrix.col(vectors[i])
      is_unique = True
      if i > 0:
        for v_u in unique_vectors:
          if v.length() < v_u.length():
            if is_approximate_integer_multiple(v, v_u):
              is_unique = False
              break
          elif is_approximate_integer_multiple(v_u, v):
            is_unique = False
            break
      if is_unique:
        unique_vectors.append(v)
        unique_indices.append(perm[i])
      i += 1
    # FIXME debugging here 
    #for direction in SST.finegrained_angles:
    #  print ('DEBUGGGG = ', direction.phi, direction.psi)

    basis_vectors = [v.elems for v in unique_vectors]
    self.candidate_basis_vectors = basis_vectors

    logger.info("Number of unique vectors: %i" % len(unique_vectors))

    for i in range(len(unique_vectors)):
      logger.debug("%s %s %s"% (
                    str(compute_functional(unique_vectors[i].elems)),
                    str(unique_vectors[i].length()),
                    str(unique_vectors[i].elems)))

    crystal_models = []
    self.candidate_basis_vectors = unique_vectors
    self.debug_show_candidate_basis_vectors()
    candidate_orientation_matrices = self.find_candidate_orientation_matrices(unique_vectors)
    crystal_model, n_indexed = self.choose_best_orientation_matrix(candidate_orientation_matrices)
    if crystal_model is not None:
      crystal_models = [crystal_model]
    else:
      crystal_models = []

    candidate_orientation_matrices = crystal_models
    self.candidate_crystal_models = candidate_orientation_matrices
