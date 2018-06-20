from psana import *

file = open("LD91_events_count","w") 

data_info = 'exp=cxid9114:dir=/reg/d/psdm/cxi/cxid9114/demo/xtc:run=95-114:smd'
ds = DataSource(data_info)

total_events = 0

for run in ds.runs():

	run_events = 0
	
    	for evt in run.events():
        	run_events = run_events + 1
	
	print 'Run: {0}: Events: {1}\n'.format(run.run(), run_events)
	file.write('Run: {0}: Events: {1}\n'.format(run.run(), run_events))
	
	total_events = total_events + run_events	

print '\nTotal events: {0}\n'.format(total_events)
file.write('\nTotal events: {0}\n'.format(total_events))

file.close()

print 'Done.'
