# Trying to be a bit smarter about the real_space_grid_search in DIALS by first doing a coarse grid search
# with the default of 0.02 grid size and then gradually reducing the grid size to 0.002 or smaller ?
# Looking to reduce the number of calls to compute functionals which is what eats up so much time

from __future__ import absolute_import, division, print_function
import abc
import math
import logging
import libtbx
from libtbx import phil

logger = logging.getLogger(__name__)

from scitbx.array_family import flex
from scitbx import matrix
from scitbx import fftpack
from cctbx import crystal, uctbx, xray

#from dials.algorithms.indexing.indexer import indexer_base, optimise_basis_vectors
#from dials.algorithms.indexing.indexer import is_approximate_integer_multiple
from dxtbx.model.experiment_list import Experiment, ExperimentList
#from dials.algorithms.indexing.real_space_grid_search import indexer_real_space_grid_search
from dials_algorithms_indexing_ext import map_centroids_to_reciprocal_space_grid
from dials.algorithms.indexing import DialsIndexError

# Import all the stuff in strategies and then write RealSpaceGridSmartSearch strategy
from dials.algorithms.indexing.basis_vector_search.strategies import Strategy, FFT1D, FFT3D, RealSpaceGridSearch
from dials.algorithms.indexing.basis_vector_search.strategies import _vector_group, _is_approximate_integer_multiple
from dials.algorithms.indexing.basis_vector_search.strategies import fft1d_phil_str, fft3d_phil_str, real_space_grid_search_phil_str


def compute_functional(vector, reciprocal_lattice_vectors):
  two_pi_S_dot_v = 2 * math.pi * reciprocal_lattice_vectors.dot(vector)
  cosines = flex.cos(two_pi_S_dot_v)
  return flex.sum(flex.cos(two_pi_S_dot_v))

real_space_grid_smart_search_phil_str = """\
coarse_sampling_grid = 0.005
  .type = float(value_min=0)
  .help = coarse sampling grid that is used to do an initial search of possible basis vectors
fine_sampling_grid = 0.0001
  .type=float(value_min=0)
  .help=fine sampling grid that further expands around the top basis vectors identified in the coarse grid search
"""


class RealSpaceGridSmartSearch(Strategy):
  """A smart real space grid search method that does an initial search based on a coarse grid and then
     expands about the top basis vectors to give more accurate values of the possible basis vectors.
  """

  phil_scope = phil.parse(real_space_grid_smart_search_phil_str)

  def __init__(self, max_cell, target_unit_cell, params=None, *args, **kwargs):
    super(RealSpaceGridSmartSearch, self).__init__(max_cell, params=params, *args, **kwargs)
    self._target_unit_cell = target_unit_cell

  def get_finegrained_SST(self, reciprocal_lattice_vectors):
    from rstbx.array_family import flex
    from rstbx.dps_core import SimpleSamplerTool
    coarse_sampling_grid=self._params.coarse_sampling_grid
    fine_sampling_grid=self._params.fine_sampling_grid
    print ('Grid sizes being used :: Coarse = %.6f Fine=%.6f'%(coarse_sampling_grid, fine_sampling_grid))
    used_in_indexing = flex.bool(reciprocal_lattice_vectors.size(), True)
    print("Indexing from %i reflections COARSE" % used_in_indexing.count(True))
    #d_min = self.params.refinement_protocol.d_min_start
    #sel = self.reflections["id"] == -1
    #if d_min is not None:
    #  sel &= 1 / self.reflections["rlp"].norms() > d_min
    #reciprocal_lattice_points = self.reflections["rlp"].select(sel)
    #print("Indexing from %i reflections COARSE" % len(reciprocal_lattice_points))


    #assert self.target_symmetry_primitive is not None
    #assert self.target_symmetry_primitive.unit_cell() is not None
    SST = SimpleSamplerTool(coarse_sampling_grid)
    SST.construct_hemisphere_grid(SST.incr)
    #cell_dimensions = self.target_symmetry_primitive.unit_cell().parameters()[:3]
    cell_dimensions = self._target_unit_cell.parameters()[:3]
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
        f = compute_functional(v.elems, reciprocal_lattice_vectors)
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
            if _is_approximate_integer_multiple(v, v_u, relative_tolerance=0.2, angular_tolerance=5.0):
              is_unique = False
              break
          elif _is_approximate_integer_multiple(v_u, v, relative_tolerance=0.2, angular_tolerance=5.0):
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

    SST.construct_hemisphere_grid_finegrained(fine_sampling_grid, coarse_sampling_grid, SST_filter)
    return SST

  def find_basis_vectors(self, reciprocal_lattice_vectors):
    """ overloading the real_space_grid_search method to use a sparser SST with highly fine grid after initially doing
        a sparse grid to evaluate grid points of interest"""


    from rstbx.array_family import flex
    from rstbx.dps_core import SimpleSamplerTool
    used_in_indexing = flex.bool(reciprocal_lattice_vectors.size(), True)
    print("Indexing from %i reflections FINE" % used_in_indexing.count(True))

    #print("Indexing from %i reflections FINE " % len(reciprocal_lattice_points))

    #Aug_Refactor
    #assert self.target_symmetry_primitive is not None
    #assert self.target_symmetry_primitive.unit_cell() is not None

    SST = self.get_finegrained_SST(reciprocal_lattice_vectors)
    #cell_dimensions = self.target_symmetry_primitive.unit_cell().parameters()[:3]
    cell_dimensions = self._target_unit_cell.parameters()[:3]
    unique_cell_dimensions = set(cell_dimensions)
    print("Number of search vectors FINE : %i"
                 %(len(SST.finegrained_angles) * len(unique_cell_dimensions)))
    vectors = flex.vec3_double()
    function_values = flex.double()
    # Do a search within the cluster of coarse grids and only return the top few values ?
    find_max_within_cluster=True
    import time
    time1=time.time()
    if find_max_within_cluster:
      top_n_values=1 # number of top scoring vectors to return in each coarse grid per unique dimension
      for count in range(len(SST.n_entries_finegrained[:-1])):
        start=SST.n_entries_finegrained[count]
        end=SST.n_entries_finegrained[count+1]
        for l in unique_cell_dimensions:
          tmp_vectors = flex.vec3_double()
          tmp_function_values = flex.double()
          for i, direction in enumerate(SST.finegrained_angles[start:end]):
            v = matrix.col(direction.dvec) * l
            f = compute_functional(v.elems, reciprocal_lattice_vectors)
            tmp_vectors.append(v.elems)
            tmp_function_values.append(f)
          perm=flex.sort_permutation(tmp_function_values, reverse=True)
          tmp_vectors=tmp_vectors.select(perm)[0:top_n_values]
          tmp_function_values=tmp_function_values.select(perm)[0:top_n_values]
          vectors.extend(tmp_vectors)
          function_values.extend(tmp_function_values)

    else:
      for i, direction in enumerate(SST.finegrained_angles):
        for l in unique_cell_dimensions:
          v = matrix.col(direction.dvec) * l
          f = compute_functional(v.elems, reciprocal_lattice_vectors)
          vectors.append(v.elems)
          function_values.append(f)
     
    time2=time.time()
    print ('FINE GRID SEARCH TIME=',time2-time1)

    perm = flex.sort_permutation(function_values, reverse=True)
    vectors = vectors.select(perm)
    function_values = function_values.select(perm)

    unique_vectors = []
    unique_indices = []
    i = 0
    time1 = time.time()
    while len(unique_vectors) < 30:
      v = matrix.col(vectors[i])
      is_unique = True
      if i > 0:
        for v_u in unique_vectors:
          if v.length() < v_u.length():
            if _is_approximate_integer_multiple(v, v_u, relative_tolerance=0.2, angular_tolerance=2.0):
              is_unique = False
              break
          elif _is_approximate_integer_multiple(v_u, v, relative_tolerance=0.2, angular_tolerance=2.0):
            is_unique = False
            break
      if is_unique:
        unique_vectors.append(v)
        unique_indices.append(perm[i])
      i += 1
    time2 = time.time()
    print ('FINE GRID UNIQUE VECTOR SEARCH TIME = ', time2-time1)
    #from IPython import embed; embed(); exit()
    # FIXME debugging here 
    #for direction in SST.finegrained_angles:
    #  print ('DEBUGGGG = ', direction.phi, direction.psi)

    basis_vectors = [v.elems for v in unique_vectors]
    self.candidate_basis_vectors = basis_vectors

    logger.info("Number of unique vectors: %i" % len(unique_vectors))

    for i in range(len(unique_vectors)):
      logger.debug("%s %s %s"% (
                    str(compute_functional(unique_vectors[i].elems, reciprocal_lattice_vectors)),
                    str(unique_vectors[i].length()),
                    str(unique_vectors[i].elems)))

    return unique_vectors, used_in_indexing
    #Aug_Refactor 
    # Stuff below is not there in strategies.py
    #crystal_models = []
    #self.candidate_basis_vectors = unique_vectors
    #self.debug_show_candidate_basis_vectors()
    #candidate_orientation_matrices = self.find_candidate_orientation_matrices(unique_vectors)
    #crystal_model, n_indexed = self.choose_best_orientation_matrix(candidate_orientation_matrices)
    #if crystal_model is not None:
    #  crystal_models = [crystal_model]
    #else:
    #  crystal_models = []
    #candidate_orientation_matrices = crystal_models
    #self.candidate_crystal_models = candidate_orientation_matrices
