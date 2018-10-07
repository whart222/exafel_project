#!/bin/bash -l

NODES=1
NUM_RANKS=$((NODES*68))

t_start=`date +%s`
RUN=24 #$1
TRIAL=6 #$2
RG=0 #$3
BASE_PATH=/global/cscratch1/sd/asmit/iota_demo/cxic0415
EXP=cxic0415

echo 'Processing LCLS run ' ${RUN}

RUN_DIR="r$(printf "%04d" ${RUN})"
TRIAL_RG="$(printf "%03d" ${TRIAL})_rg$(printf "%03d" ${RG})"
RG_3d="rg_$(printf "%03d" ${RG})"

mkdir -p ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/out/
mkdir -p ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/stdout/

cp ${BASE_PATH}/input/process_batch.phil ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/params_1.phil

#export PMI_MMAP_SYNC_WAIT_TIME=600
srun -n ${NUM_RANKS} -c 4 --cpu_bind=cores shifter ${BASE_PATH}/docker_xtc_process.sh ${RUN} ${TRIAL} ${RG} ${BASE_PATH} ${EXP} production 0

#srun -n 1 -c 68 --cpu_bind=cores shifter ${BASE_PATH}/docker_single_ts.sh ${RUN} ${TRIAL} ${RG} ${BASE_PATH} ${EXP} production 0

#mkdir ${PWD}/tmp
#cp  ./input/process_batch.phil ${PWD}/tmp/process_batch.phil
#cp  ./xtc_process.py ${PWD}/tmp/xtc_process.py
#srun -n 1088 -c 4 --cpu_bind=cores shifter ./index_single.sh cxic0415 24 0 debug 0 ${PWD}/tmp
t_end=`date +%s`
echo IOTA_XTC_JobCompleted_TimeElapsed $((t_end-t_start)) $t_start $t_end ${RUN}
