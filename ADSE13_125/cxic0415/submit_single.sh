#!/bin/bash -l
t_start=`date +%s`

#export PMI_MMAP_SYNC_WAIT_TIME=600
#mkdir -p /global/cscratch1/sd/asmit/iota_demo/cxic0415/${RUN_DIR}/${TRIAL_RG}/out
#mkdir -p /global/cscratch1/sd/asmit/iota_demo/cxic0415/${RUN_DIR}/${TRIAL_RG}/stdout
cp  ./input/process_batch.phil ${PWD}/tmp/process_batch.phil
cp  ./xtc_process.py ${PWD}/tmp/xtc_process.py
srun -n 8 -c 4 --cpu_bind=cores shifter ./index_single.sh cxic0415 1 0 debug 4 ${PWD}/tmp
t_end=`date +%s`
echo PSJobCompleted TotalElapsed $((t_end-t_start)) $t_start $t_end
