from __future__ import division
import sys
from libtbx import easy_pickle
from scitbx.matrix import sqr, col
from cctbx import crystal # dependency for integration pickles
from cctbx_orientation_ext import crystal_orientation

"""
Script that examines a set of cctbx.xfel experiments and writes out their basis vectors
in reciprocal space and/or realspace  in gnuplot format.

Usage: Supply pickle file with a list of experiments on the command line. Create arrows.p file. Run gnuplot, then enter load "arrows.p".
"""

f = open("arrows.p",'w')

init_str = '''set xrange [-50:50]
set yrange [-50:50]
set zrange [-50:50]
set xlabel "x*" offset 0
set ylabel "y*" offset 0
set zlabel "z*" offset 0
#set xtics offset -5
#set ytics offset 5
set format x ""
set format y ""
set format z ""
set ticslevel 1.0
set border lw 4
splot sqrt(-1) \n'''
f.write(init_str)
filename = sys.argv[1]
try:
  data = easy_pickle.load(filename)
except Exception as e:
  print "Couldn't read", filename
is_reciprocal = True
#from IPython import embed; embed(); exit()
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
    #from IPython import embed; embed(); exit()
 #   assert  ori.unit_cell().parameters()==ori_best.unit_cell().parameters(), 'unit cell '
  A = sqr(ori_best.reciprocal_matrix()).transpose().inverse()
  abasis = A * col((1,0,0))
  bbasis = A * col((0,1,0))
  cbasis = A * col((0,0,1))
  print 'ANGLES A-B=%5.3f, B-C=%5.3f, C-A=%5.3f' %(abasis.angle(bbasis,deg=True), bbasis.angle(cbasis,deg=True), cbasis.angle(abasis,deg=True))
  print A[0]
    #from IPython import embed; embed(); exit()
  f.write("set arrow to % 6.6f, % 6.6f, % 6.6f lc rgb 'red'  \n"%(abasis[0],abasis[1],abasis[2]))
  f.write("set arrow to % 6.6f, % 6.6f, % 6.6f lc rgb 'green'\n"%(bbasis[0],bbasis[1],bbasis[2]))
  f.write("set arrow to % 6.6f, % 6.6f, % 6.6f lc rgb 'blue' \n"%(cbasis[0],cbasis[1],cbasis[2]))
f.close()
