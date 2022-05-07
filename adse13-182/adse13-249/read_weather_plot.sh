#! /bin/bash

# example output from weather2.py:

#OK
#rank_3.log 1
#the median weather time is 1.35002 for job 2122552
#Five number summary of 32 good image processing times: ['1.33821', '1.34492', '1.34991', '1.36550', '1.77553']
#The total envelope time is 11.3 seconds
#The total MPI communicator time is 29.9 seconds, with 18.7 sec before 'foreach' and -0.1 sec trailing
#The total Python time is 63.4 seconds, with 33.5 sec for imports and 0.0 sec trailing
#B: Python imports  33.50
#C: MPI gather SF   14.37
#D: MPI broadcast    3.96
#E: logger redirect  0.00, mean   0.00
#F: set CUDA device  0.34
#G: big data to GPU  0.37, mean   0.37

export envelope_time=`grep "envelope" weather.out | awk '{ print $6 }'`
export python_time=`grep "total Python" weather.out | awk '{ print $6 }'`
export comm_time=`grep "MPI communicator" weather.out | awk '{ print $7 }'`
export weather_time=`grep "weather time" weather.out | awk '{ print $6 }'`

echo env$'\t'py$'\t'comm$'\t'weather
echo $envelope_time,$'\t'$python_time,$'\t'$comm_time,$'\t'$weather_time
