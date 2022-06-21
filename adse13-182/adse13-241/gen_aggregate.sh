#! /bin/bash

echo "env, py, comm, weather, jobid, nodes, ranks, r/gpu, evts" > aggregate.out
for infostr in $(cat record)
  do #echo $infostr
	#info
  jobid=${infostr%_N*}
  wfile="${jobid}/weather.out"
  tmp=${infostr%_gpn*}
  nodes=${tmp#*_N}
  tmp=${infostr%_n*_nimg*}
  gpn=${tmp#*_gpn}
  tmp=${infostr%_nimg*}
  ranks=${tmp#*_n}
  evts=${infostr#*_nimg}
  rpgpu=$((ranks / (nodes * gpn)))
  #echo $jobid nodes: $nodes gpn: $gpn ranks: $ranks rpgpu: $rpgpu evts: $evts
  #timings
  env=`grep 'envelope time' $wfile | awk '{ print $6 }'`
	py=`grep 'total Python time' $wfile | awk '{ print $6 }'`
	comm=`grep 'communicator time' $wfile | awk '{ print $7 }'`
	weather=`grep 'weather time' $wfile | awk '{ print $6 }'`
	echo "$env, $py, $comm, $weather, $jobid, $nodes, $ranks, $rpgpu, $evts" >> aggregate.out
done

#2363201 N1_gpn4_n4_nimg32
#2363225 N1_gpn4_n4_nimg320


#env     py      comm    weather         jobid           nodes   ranks   r/gpu   evts
#11.3,   76.7,   33.9,   1.36154,        2122549,        1,      4,      1,      32


#OK
#rank_25.log 1
#the median weather time is 1.23379 for job 2363438
#Five number summary of 3200 good image processing times: ['1.05761', '1.23137', '1.23379', '1.24171', '2.06225']
#The total envelope time is 53.6 seconds
#The total MPI communicator time is 64.2 seconds, with 11.2 sec before 'foreach' and -0.6 sec trailing
#The total Python time is 74.4 seconds, with 10.2 sec for imports and 0.0 sec trailing
#B: Python imports  10.23
#C: MPI gather SF    4.84
#D: MPI broadcast    5.49
#E: logger redirect  0.01, mean   0.01
#F: set CUDA device  0.88
#G: big data to GPU  0.22, mean   0.19
