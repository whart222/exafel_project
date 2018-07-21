from __future__ import division
import sys, os
from libtbx import easy_pickle
from scitbx.matrix import sqr, col
from cctbx import crystal # dependency for integration pickles
from scitbx.math import flex
from cctbx_orientation_ext import crystal_orientation

"""
script to print out the angle in degrees between basis vectors of a crystal
Input is an experimental model(.json)
"""

filename = sys.argv[1]
is_reciprocal = True
from dxtbx.model.experiment_list import ExperimentListFactory
exp = ExperimentListFactory.from_json_file(filename, check_format=False)

ori = crystal_orientation(exp.crystals()[0].get_A(), is_reciprocal)
A = sqr(ori.reciprocal_matrix())
abasis = A * col((1,0,0))
bbasis = A * col((0,1,0))
cbasis = A * col((0,0,1))
print 'ANGLES A-B=%5.3f, B-C=%5.3f, C-A=%5.3f' %(abasis.angle(bbasis,deg=True), bbasis.angle(cbasis,deg=True), cbasis.angle(abasis,deg=True))
print 'A-matrix = ',A
