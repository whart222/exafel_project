#!/bin/bash -l

NODES=32
NUM_RANKS=$((NODES*68))


# Process run group 37 from LS49
for i in `seq 223 223`
do
  echo 'Processing LCLS run ' $i
  RUN=$i
  TRIAL=2
  RG=37
  BASE_PATH=/global/cscratch1/sd/asmit/LS49/iota_v2_runs_165_255
  EXP=mfxls4916

  RUN_DIR="r$(printf "%04d" ${RUN})"
  TRIAL_RG="$(printf "%03d" ${TRIAL})_rg$(printf "%03d" ${RG})"
  RG_3d="rg_$(printf "%03d" ${RG})"

  #rm -r ${BASE_PATH}/${RUN_DIR}
  mkdir -p ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/out
  mkdir -p ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/stdout

  cp ${BASE_PATH}/phil_files/${RG_3d}.phil ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/params_1.phil


  #sbcast -f ${BASE_PATH}/phil_files/${RG_3d}.phil /tmp/params_1.phil
  #sbcast -f /global/homes/a/asmit/cctbx.xfel/modules/exafel_project/ADSE13_25/xtc_process_iota_baseline.py /tmp/xtc_process_iota_baseline.py
  #sbcast -f /global/cscratch1/sd/asmit/LS49/iota/masks/mask.pickle /tmp/mask.pickle
  #sbcast -p /global/homes/a/asmit/cctbx.xfel/modules/cctbx_project/xfel/command_line/xtc_process.py /tmp/xtc_process.py


  t_start=`date +%s`
  srun -n ${NUM_RANKS} -c 4 --cpu_bind=cores shifter ${BASE_PATH}/docker_xtc_process.sh ${RUN} ${TRIAL} ${RG} ${BASE_PATH} ${EXP}
  #srun -n ${NUM_RANKS} -c 4 --cpu_bind=cores shifter ${BASE_PATH}/mpi_test.sh ${RUN} ${TRIAL} ${RG} ${BASE_PATH} ${EXP}
  t_end=`date +%s`
  echo IOTA_XTC_JobCompleted_TimeElapsed $((t_end-t_start)) $t_start $t_end ${RUN}


done
