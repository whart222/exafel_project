from __future__ import absolute_import,print_function, division
import matplotlib.pyplot as plt
import sys,os
from iotbx.detectors.cspad_detector_formats import reverse_timestamp
from libtbx.phil import parse
from libtbx.utils import Sorry

message = ''' script to get a sense of the computational performance of every rank while processing xtc streams
              End product is a plot of wall time vs MPI rank number with every data point being that of a frame
              processed by xtc_process. The information is read in from the debug files created by xtc_process
              Example usage on cxic0415 processed demo data -
              libtbx.python analyze_computational_performance.py
              input_path=/global/project/projectdirs/lcls/mona/output_demo18/cxic0415/output/debug
'''
phil_scope = parse('''
  input_path = .
    .type = str
    .help = path to where the debug txt files are located
  num_nodes = 100
    .type = int
    .help = Number of nodes used to do data processing. Used in timing information
  num_cores_per_node = 68
    .type = int
    .help = Number of cores per node in the machine (default is for Cori KNL)
  wall_time = 3600
    .type = int
    .help = total wall time (seconds) taken for job to finish. Used for plotting node-partitioning
  plot_title = Computational weather plot
    .type = str
    .help = title of the computational weather plot
  pickle_plot = False
    .type = bool
    .help = If True, will pickle matplotlib session so that it can be opened later for analysis/viewing \
            https://stackoverflow.com/questions/29160177/matplotlib-save-file-to-be-reedited-later-in-ipython
  pickle_filename = fig_object.pickle
    .type = str
    .help = Default name of pickled matplotlib plot saved to disk
''')

def params_from_phil(args):
  user_phil = []
  for arg in args:
    if os.path.isfile(arg):
      user_phil.append(parse(file_name=arg))
    else:
      try:
        user_phil.append(parse(arg))
      except Exception as e:
        raise Sorry("Unrecognized argument: %s"%arg)
  params = phil_scope.fetch(sources=user_phil).extract()
  return params

def run(params):
  counter = 0
  reference = None
  root=params.input_path
  fig_object = plt.figure()
  for filename in os.listdir(root):
    if os.path.splitext(filename)[1] != '.txt': continue
    if 'debug' not in filename: continue
    timepoints = []
    rank = int(filename.split('_')[1].split('.')[0])
    counter += 1
    print (filename)
    for line in open(os.path.join(root,filename)):
      hostname, psanats, ts, status, result = line.strip().split(',')
      if reference is None:
        sec, ms = reverse_timestamp(ts)
        reference = sec+ms*1e-3

      if status in ['stop','done']:
        sec, ms = reverse_timestamp(ts)
        timepoints.append((sec + ms*1.e-3)-reference)
        ok = True
      else:
        ok = False
    plt.plot(timepoints, [rank]*len(timepoints), 'b.')
    if not ok:
      sec, ms = reverse_timestamp(ts)
      plt.plot([(sec+ms*1e-3) - reference], [rank], 'rx')
    #if counter > 100: break
  for i in range(params.num_nodes):
    plt.plot([0,params.wall_time], [i*params.num_cores_per_node-0.5, i*params.num_cores_per_node-0.5], 'r-')
  plt.xlabel('Wall time (sec)')
  plt.ylabel('MPI Rank Number')
  plt.title(params.plot_title)
  if params.pickle_plot:
    from libtbx.easy_pickle import dump
    dump('%s'%params.pickle_filename, fig_object)
  plt.show()

if __name__ == '__main__':
  if '--help' in sys.argv[1:] or '-h' in sys.argv[1:]:
    print (message)
    exit()
  params = params_from_phil(sys.argv[1:])
  run(params)
