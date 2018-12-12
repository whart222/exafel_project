from __future__ import division, print_function, absolute_import

message = '''Pass in 2 half dataset mtz file and it should return CC_anom value'''


from iotbx import reflection_file_reader
from cctbx.array_family import flex
import sys

def run(mtz1, mtz2):
  # read in mtz files
  reflection_file_1 = reflection_file_reader.any_reflection_file(mtz1)
  miller_arrays_1 = reflection_file_1.as_miller_arrays()
  reflection_file_2 = reflection_file_reader.any_reflection_file(mtz2)
  miller_arrays_2 = reflection_file_2.as_miller_arrays()

  ma_anom_diff_1 = miller_arrays_1[0].anomalous_differences()
  ma_anom_diff_2 = miller_arrays_2[0].anomalous_differences()
  ma_anom_diff_1.show_summary()
  ma_anom_diff_2.show_summary()
  ma_anom_diff_1_cs = ma_anom_diff_1.common_set(ma_anom_diff_2)
  ma_anom_diff_2_cs = ma_anom_diff_2.common_set(ma_anom_diff_1)
  # Make sure the 2 arrays have the same size
  cc_anom = flex.linear_correlation(ma_anom_diff_1_cs.data(), ma_anom_diff_2_cs.data()).coefficient()
  print ('Value of CC_anom for the dataset is = ',cc_anom)
  #from IPython import embed; embed(); exit()

if __name__ == '__main__':
  print (message)
  run(sys.argv[1], sys.argv[2])
