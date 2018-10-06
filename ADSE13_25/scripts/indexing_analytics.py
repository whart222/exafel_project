from __future__ import print_function
import sys, os, math
from iotbx.detectors.cspad_detector_formats import reverse_timestamp
from matplotlib import pyplot as plt
from dials.array_family import flex

''' 
Code for calculating analytics for indexing

Run code like libtbx.python indexing_analytics.py {input_path} {number_of_nodes} {out_logfile}
input_path  : Path where pickle files & debug files exist

input_path should be like that used at LCLS i.e full path to {run_number}/{trial_rg}
Assumes folder structure


Input Path    -----> out -----> debug folder
                         -----> pickle and json files
              -----> stdout


number_of_nodes = Number of nodes used to do data processing
 		  Used in timing information


out_logfile = (optional). Needed for picking up timing information

'''



num_nodes = int(sys.argv[2])
root = os.path.join(sys.argv[1], 'out')
run_num = '99999'
steps_d = {}
debug_root = os.path.join(root, 'debug')
out_logfile = None
if len(sys.argv) == 4:
  out_logfile = os.path.join(sys.argv[1], 'stdout', sys.argv[3])

debug_mode=False

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

event_number = 0
hits = []
events_list = []
indexing_time = {}
indexing_time_all = {} # Both for failed and successful indexing trials
for filename in os.listdir(debug_root):
  if os.path.splitext(filename)[1] != ".txt": continue
#  print filename
  current_ts = ""
  current_reverse_ts = ""
  prev_step = ""
  prev_time = None
  with open(os.path.join(debug_root, filename)) as logfile:
    for line in logfile:
      try:
        hostname, ts, now, status, step = line.split(',')
      except ValueError:
        print (line); raise
      now_s, now_ms = reverse_timestamp(now)
      now = now_s + (1e-3 * now_ms)
      if step.strip() == 'start':
        event_number += 1
        curr_ts, curr_ms = reverse_timestamp(ts)
        current_reverse_ts = curr_ts + (1e-3*curr_ms)
        #if '2018-05-01T14:50Z21.976' == ts:
        #  print (line + 'THIS MIGHT BE A DUPLICATE\n')
        events_list.append(ts)
        if prev_time is not None:
          #print prev_step.strip(), "took", now - prev_time, "seconds"
	  add_step(prev_step, now-prev_time)
      else:
        #print prev_step.strip(), "took", now - prev_time, "seconds"
	add_step(prev_step, now-prev_time)
     

      if prev_step.strip() == 'index_start':
        indexing_time_all[ts] = now-prev_time
       
        #print ('index_start recognized and recorded')

      if step.strip() == 'index_start':
        hits.append(ts) 
      if step.strip() == 'refine_start':
        indexing_time[ts]=now-prev_time
      
      current_ts = ts
      prev_step = step
      prev_time = now

#print steps_d.keys()
#print indexing_time.keys()
#print len(hits)
total_idx_time=0
for event in indexing_time.keys():
#  print event, indexing_time[event]
  total_idx_time +=indexing_time[event]

recorded_hits = []
idx_attempt_time = []
idx_successful_time = []

#fout = open('indexing_timing_' + run_num+'.dat','w')
#fout.write('Event Number             hits       indexed         t_indexed              t_indexed_attempted  \n' )
for ii,event in enumerate(events_list):
  if event in hits:
    is_hit=1
    if event not in recorded_hits:
      recorded_hits.append(event)
    else:
      if debug_mode:
        print ('Duplicate Event ? = ', event)
  else:
    is_hit=0
  if event in indexing_time:
    is_idx=1
    t_idx=indexing_time[event]
    idx_successful_time.append(t_idx)
    assert event in indexing_time_all, 'Event not present in indexing_time_all'
    t_idx2 = indexing_time_all[event]
    idx_attempt_time.append(t_idx2)
  elif event in indexing_time_all:
    is_idx=0.5
    assert event not in indexing_time, 'Event should not be present in indexing_time'
    t_idx=0.0
    t_idx2 = indexing_time_all[event]
    idx_attempt_time.append(t_idx2)
  else:
    is_idx=0
    t_idx=0.0
    t_idx2 =0.0
  
#  fout.write('%s  %3.1f  %3.1f  %12.7f  %12.7f\n' %(event, is_hit, is_idx, t_idx, t_idx2)) 
#fout.close()

# Extract timing information from log file
if out_logfile is not None:
  total_time = []
  run_number = int(root.strip().split('/')[-3][1:])
  print (run_number)
  with open(out_logfile, 'r') as flog:
    for line in flog:
      if 'IOTA_XTC_SingleRank_TimeElapsed' in line:
        ax = line.split()
        if int(ax[-1]) == run_number:
          total_time.append(float(ax[1]))


  node_hours = max(total_time)*num_nodes/3600.0


# Unit cell and RMSD statistics for that run
all_uc_a = flex.double()
all_uc_b = flex.double()
all_uc_c = flex.double()
all_uc_alpha = flex.double()
all_uc_beta = flex.double()
all_uc_gamma = flex.double()

from dxtbx.model.experiment_list import ExperimentListFactory
from dials.algorithms.refinement.prediction import ExperimentsPredictor
from libtbx.easy_pickle import load
from scitbx.matrix import col
dR = flex.double()
for filename in os.listdir(root):
  if 'integrated_experiments' not in os.path.splitext(filename)[0] or os.path.splitext(filename)[1] != ".json": continue

  fjson=os.path.join(root, filename)
  experiments = ExperimentListFactory.from_json_file(fjson)
  for crystal in experiments.crystals():
    all_uc_a.append(crystal.get_unit_cell().parameters()[0]) 
    all_uc_b.append(crystal.get_unit_cell().parameters()[1]) 
    all_uc_c.append(crystal.get_unit_cell().parameters()[2]) 
    all_uc_alpha.append(crystal.get_unit_cell().parameters()[3]) 
    all_uc_beta.append(crystal.get_unit_cell().parameters()[4]) 
    all_uc_gamma.append(crystal.get_unit_cell().parameters()[5]) 

  fpickle = os.path.join(root, filename.split('integrated_experiments')[0]+'indexed.pickle')
  reflections = load(fpickle)
  ref_predictor = ExperimentsPredictor(experiments, force_stills=experiments.all_stills())
  reflections = ref_predictor(reflections)
  for refl in reflections:
    dR.append((col(refl['xyzcal.mm']) - col(refl['xyzobs.mm.value'])).length())

# Now print out all relevant statistics
print ('-'*80)
print ('|'+' '*80+'|\n'+'|'+ ' '*20 + 'Analytics Package for Indexing'+' '*30+'|\n|'+' '*80+'|')
print ('-'*80)
print ('Getting stats for data in : ',root)
print ('====================== Indexing and Timing Statistics ============================')
print ('Number of Hits = ', len(hits))
print ('Number of images successfully indexed = ', len(idx_successful_time))
print ('Total time spent in indexing (hrs) = ', sum(idx_attempt_time)/3600.0)
print ('Time spent in indexing successfully (hrs) = ', sum(idx_successful_time)/3600.0)
if out_logfile is not None:
  print ('Total Node-hours with %d nodes = %.2f (hrs)'%(num_nodes, node_hours))
print ('====================== Unit Cell & RMSD Statistics ============================')
print ('a-edge (A) : %.2f +/- %.2f' % (flex.mean(all_uc_a),flex.mean_and_variance(all_uc_a).unweighted_sample_standard_deviation()))
print ('b-edge (A) : %.2f +/- %.2f' % (flex.mean(all_uc_b),flex.mean_and_variance(all_uc_b).unweighted_sample_standard_deviation()))
print ('c-edge (A) : %.2f +/- %.2f' % (flex.mean(all_uc_c),flex.mean_and_variance(all_uc_c).unweighted_sample_standard_deviation()))
print ('alpha (deg) : %.2f +/- %.2f' % (flex.mean(all_uc_alpha),flex.mean_and_variance(all_uc_alpha).unweighted_sample_standard_deviation()))
print ('beta (deg) : %.2f +/- %.2f' % (flex.mean(all_uc_beta),flex.mean_and_variance(all_uc_beta).unweighted_sample_standard_deviation()))
print ('gamma (deg) : %.2f +/- %.2f' % (flex.mean(all_uc_gamma),flex.mean_and_variance(all_uc_gamma).unweighted_sample_standard_deviation()))
print ('Total RMSD i.e calc - obs for Bragg spots (um) = ', 1000.0*math.sqrt(dR.dot(dR)/len(dR)))
#from IPython import embed; embed(); exit()


