import sys, os
from iotbx.detectors.cspad_detector_formats import reverse_timestamp
import numpy as np
from matplotlib import pyplot as plt

# Run this script in the debug directory as 
# libtbx.python histogram_timings.py .       #NOTE THE DOT !!!



root = os.getcwd() #sys.argv[1]


def add_step(step, duration, steps_d,ts,fname,now_step):
  step = step.strip()
  my_step = step
  if step.endswith("_start"):
    step = "_".join(step.split('_')[:-1])
  if any([s in step for s in ['_ok_','_failed_']]):
    step = "_".join(step.split('_')[:-1])
  if 'not_enough_spots' in step:
    step = "_".join(step.split('_')[:-1])

  if step not in steps_d:
    steps_d[step] = []
  steps_d[step].append(duration)
  if 'indexing_failed' in my_step:
    #print 'psanagpu999,%s,%s,fail'%(ts,ts)
    #print 'psanagpu999',ts,ts,'fail'
    print 'F_D_B',ts,' ',duration,fname,my_step.strip()
  return steps_d
  

def get_timing_info(root, rank=0, size=1):
  steps_d = {}
  for ii,filename in enumerate(os.listdir(root)):
    if os.path.splitext(filename)[1] != ".txt": continue
    if (ii+rank)%size !=0: continue
    current_ts = ""
    prev_step = ""
    prev_time = None
    with open(os.path.join(root, filename)) as logfile:
      for line in logfile:
        try:
          hostname, ts, now, status, step = line.split(',')
        except ValueError:
          print 'ERROR',line; raise
        now_s, now_ms = reverse_timestamp(now)
        now = now_s + (1e-3 * now_ms)
        #from IPython import embed; embed(); exit()
        if step.strip() == 'start':
          if prev_time is not None:
            #print prev_step.strip(), "took", now - prev_time, "seconds"
  	    steps_d=add_step(prev_step, now-prev_time,steps_d,ts,filename,step)
        else:
          #print prev_step.strip(), "took", now - prev_time, "seconds"
  	  steps_d=add_step(prev_step, now-prev_time,steps_d,ts,filename,step)
        current_ts = ts
        prev_step = step
        prev_time = now
    #if 'index' in step:
    #  print 'FWWWWW',step,filename,ts
  
  #print steps_d.keys()
  return steps_d

def run(root):
  comm=None
  rank=0
  if True:
    try:
      from mpi4py import MPI
    except ImportError:
      raise Sorry("MPI not found")
    comm=MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    #print (rank, size)

    steps_d_mpi=get_timing_info(root,rank=rank,size=size) 
    steps_d_mpi=comm.gather(steps_d_mpi, root=0)
    if rank !=0: return
    #from IPython import embed; embed(); exit()
    steps_d = {}
    for entry in steps_d_mpi:
      for key in entry:
        if key not in steps_d:
          steps_d[key] = []
        steps_d[key].extend(entry[key])

  exit()
  for key in steps_d:
    plt.figure()
    plt.title(key)
    plt.hist(steps_d[key], bins=100)
  plt.show()

if __name__=='__main__':
  run(root)
