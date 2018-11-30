import sys, os
from iotbx.detectors.cspad_detector_formats import reverse_timestamp
import numpy as np
from matplotlib import pyplot as plt

# Run this script in the debug directory as 
# libtbx.python histogram_timings.py .       #NOTE THE DOT !!!



root = os.getcwd() #sys.argv[1]

steps_d = {}

def add_step(step, duration):
  step = step.strip()
  if step.endswith("_start"):
    step = "_".join(step.split('_')[:-1])
  if any([s in step for s in ['_ok_','_failed_']]):
    step = "_".join(step.split('_')[:-1])
  if 'not_enough_spots' in step:
    step = "_".join(step.split('_')[:-1])

  if step not in steps_d:
    steps_d[step] = []
  steps_d[step].append(duration)

for filename in os.listdir(root):
  if os.path.splitext(filename)[1] != ".txt": continue
  print filename
  current_ts = ""
  prev_step = ""
  prev_time = None
  with open(os.path.join(root, filename)) as logfile:
    for line in logfile:
      try:
        hostname, ts, now, status, step = line.split(',')
      except ValueError:
        print line; raise
      now_s, now_ms = reverse_timestamp(now)
      now = now_s + (1e-3 * now_ms)
      from IPython import embed; embed(); exit()
      if step.strip() == 'start':
        if prev_time is not None:
          print prev_step.strip(), "took", now - prev_time, "seconds"
	  add_step(prev_step, now-prev_time)
      else:
        print prev_step.strip(), "took", now - prev_time, "seconds"
	add_step(prev_step, now-prev_time)
      current_ts = ts
      prev_step = step
      prev_time = now

print steps_d.keys()

for key in steps_d:
  plt.figure()
  plt.title(key)
  plt.hist(steps_d[key], bins=100)
plt.show()
