#!/bin/bash

export SIT_PSDM_DATA=/global/cscratch1/sd/asmit
export SIT_DATA=/global/cscratch1/sd/asmit/LS10/psdm/psdm/data
source /build/setpaths.sh
#
RUN=$1
TRIAL=$2
RG=$3
BASE_PATH=$4
EXP=$5

RUN_DIR="r$(printf "%04d" ${RUN})"
TRIAL_RG="$(printf "%03d" ${TRIAL})_rg$(printf "%03d" ${RG})"
RG_3d="rg_$(printf "%03d" ${RG})"

START_XTC=$(date +"%s")

#time libtbx.python /modules/exafel_project/ADSE13_25/xtc_process_iota_baseline.py output.output_dir=${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/out output.logging_dir=None input.trial=${TRIAL} input.experiment=${EXP} input.run_num=${RUN} input.cfg=None input.xtc_dir=None input.use_ffb=False input.rungroup=${RG} mp.method=mpi /tmp/params_1.phil
#time libtbx.python /tmp/xtc_process_iota_baseline.py output.output_dir=${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/out output.logging_dir=None input.trial=${TRIAL} input.experiment=${EXP} input.run_num=${RUN} input.cfg=None input.xtc_dir=None input.use_ffb=False input.rungroup=${RG} mp.method=mpi /tmp/params_1.phil

time libtbx.python /modules/exafel_project/ADSE13_25/xtc_process_iota_baseline.py output.output_dir=${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/out output.logging_dir=None input.trial=${TRIAL} input.experiment=${EXP} input.run_num=${RUN} input.cfg=None input.xtc_dir=None input.use_ffb=False input.rungroup=${RG} mp.method=mpi ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/params_1.phil
#time libtbx.python /global/homes/a/asmit/cctbx.xfel/modules/cctbx_project/xfel/command_line/xtc_process.py output.output_dir=${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/out output.logging_dir=None input.trial=${TRIAL} input.experiment=${EXP} input.run_num=${RUN} input.cfg=None input.xtc_dir=None input.use_ffb=False input.rungroup=${RG} mp.method=mpi dispatch.max_events=100 ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/params_1.phil

END_XTC=$(date +"%s")
ELAPSED=$((END_XTC-START_XTC))
echo IOTA_XTC_SingleRank_TimeElapsed ${ELAPSED} ${START_XTC} ${END_XTC} ${RUN}
