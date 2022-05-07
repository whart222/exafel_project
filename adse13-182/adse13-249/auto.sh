#! /bin/bash
# execute me from any of the $PSCRATCH/nesap/adse13-249/work/jobid directories
# to begin making weather plots and printing timing information for each
# requested jobid supplied on the command line.

while [ 1 ]
do
read -p 'jobid: ' jobid
cd ../$jobid
libtbx.python ../../current_scripts/weather2.py > weather.out
./../../current_scripts/read_weather_plot.sh
done
