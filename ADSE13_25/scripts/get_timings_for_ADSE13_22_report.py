'''
Get timings to be reported in ADSE13-22 IOTA report for DOE
'''
from __future__ import division


import sys
num_nodes=64
time_dict = {}
with open(sys.argv[1]+'.log','r') as fin:
  for line in fin:
    if 'IOTA_XTC_SingleRank_TimeElapsed' in line:
      ax = line.split()
      time_secs=float(ax[1])
      time_hrs=time_secs/3600.0
      run_num = int(ax[-1])
      if run_num not in time_dict:
        time_dict[run_num] = []
      time_dict[run_num].append(time_hrs*num_nodes)


with open('timings_srs.dat','w') as fout:
  for entry in time_dict:
    fout.write('%6d  %8.3f\n'%(entry, max(time_dict[entry])))
