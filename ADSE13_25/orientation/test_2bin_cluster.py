from __future__ import division
import sys
from libtbx import easy_pickle
from scitbx.matrix import sqr, col
from cctbx import crystal # dependency for integration pickles
from cctbx_orientation_ext import crystal_orientation

"""
Script to analysze the A matrices from multiple experimental models saved in a pickle file. This includes clustering
"""
from exafel_project.ADSE13_25.consensus_functions import get_uc_consensus as get_consensus
filename = sys.argv[1]
try:
  data = easy_pickle.load(filename)
except Exception, e:
  print "Couldn't read", filename

known_crystal_models = get_consensus(data, show_plot=True)

#from IPython import embed; embed(); exit()
is_reciprocal = True
ori0 = crystal_orientation(data[0].crystals()[0].get_A(), is_reciprocal)

for i,exp in enumerate(data):
  ori = crystal_orientation(exp.crystals()[0].get_A(), is_reciprocal)
  try:
    best_similarity_transform = ori.best_similarity_transformation(
      other = ori0, fractional_length_tolerance = 1.00,
      unimodular_generator_range=1)
    ori_best=ori.change_basis(best_similarity_transform)
  except Exception:
    ori_best = ori
  A = sqr(ori_best.reciprocal_matrix()).transpose().inverse()
  abasis = A * col((1,0,0))
  bbasis = A * col((0,1,0))
  cbasis = A * col((0,0,1))
  print A[0]
