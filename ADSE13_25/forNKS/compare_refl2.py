from __future__ import division
from libtbx.easy_pickle import load
from scitbx.matrix import col
from scitbx.math import five_number_summary

message = ''' this script compares predicted (x,y) vs observed (x,y) on the detector '''
print (message)

def apply_filter(hkl_tuple, filter_array = [1,1,1]):
  return tuple((hkl_tuple[0]*filter_array[0],hkl_tuple[1]*filter_array[1] , hkl_tuple[2]*filter_array[2]))


refl_iota = load('idx-step5_MPIbatch_000064.img_indexed.pickle')

iota_dr = []

for ii in range(len(refl_iota)):
  xyzobs_iota = refl_iota['xyzobs.px.value'][ii]
  xyzcal_iota = refl_iota['xyzcal.px'][ii]

  iota_dr.append((col(xyzobs_iota)-col(xyzcal_iota)).length())

print ('Now analyzing: Printing 5-number summary of dR = |robs-rcal|')
print (five_number_summary(iota_dr))
print ('Now plotting histogram of difference in dR = |robs-rcal|')
import matplotlib.pyplot as plt
plt.figure(1)
plt.hist(iota_dr,bins=20)
#plt.xlim([-1, max(max(base_dr), max(iota_dr))])
plt.show()


#from IPython import embed; embed(); exit()
