#!/bin/bash

source /cctbx/build/setpaths.sh
ts="2015-07-26T12:43Z45.065"

RUN=$1
TRIAL=$2
RG=$3
BASE_PATH=$4
EXP=$5
CMDMODE=$6
LIMIT=$7

RUN_DIR="r$(printf "%04d" ${RUN})"
TRIAL_RG="$(printf "%03d" ${TRIAL})_rg$(printf "%03d" ${RG})"
RG_3d="rg_$(printf "%03d" ${RG})"

START_XTC=$(date +"%s")

# mask and metrology files are from the current dir
IN_DIR=${BASE_PATH}/input
#DATA_DIR=${BASE_PATH}/../small_xtc2/${EXP}
DATA_DIR=/global/cscratch1/sd/asmit/iota_demo/small_xtc2/${EXP}

export PS_CALIB_DIR=$IN_DIR
export PS_SMD_N_EVENTS=1
export PS_SMD_NODES=10

cctbx_args="input.experiment=${EXP} input.run_num=${RUN} output.logging_dir=None output.output_dir=${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/out format.cbf.invalid_pixel_mask=${IN_DIR}/mask.pickle ${BASE_PATH}/${RUN_DIR}/${TRIAL_RG}/params_1.phil dump_indexed=True input.xtc_dir=${DATA_DIR} "

if [ "${CMDMODE}" = "debug_timestamp" ]; then
  cctbx_args="$cctbx_args debug.event_timestamp=${ts}"
fi

#echo ${cctbx_args}
#exit


if [ "$LIMIT" -ne 0 ]; then
    cctbx_args="$cctbx_args max_events=${LIMIT}"
fi

if [ "${CMDMODE}" = "pythonprof" ]; then
    libtbx.python -m cProfile -s tottime libtbx.python /cctbx/modules/exafel_project/ADSE13_25/xtc_process_iota_srs.py ${cctbx_args}
  
elif [ "${CMDMODE}" = "strace" ]; then
    strace -ttt -f -o $$.log libtbx.python /cctbx/modules/exafel_project/ADSE13_25/xtc_process_iota_srs.py ${cctbx_args}

elif [ "${CMDMODE}" = "debug" ]; then
    #libtbx.python ${PWD}/xtc_process_iota_srs_ps2.py ${cctbx_args}
    libtbx.python ${PWD}/xtc_process.py ${cctbx_args}

elif [ "${CMDMODE}" = "debug_timestamp" ]; then
    #python ${PWD}/xtc_process_iota_srs_ps2.py ${cctbx_args}
    libtbx.python ${PWD}/xtc_process.py ${cctbx_args}
  
else
    echo "Running in production mode"
    #libtbx.python /cctbx/modules/exafel_project/ADSE13_25/xtc_process_iota_srs.py ${cctbx_args}
    #libtbx.python /cctbx/modules/cctbx_project/xfel/command_line/xtc_process.py ${cctbx_args}
    cctbx.xfel.xtc_process ${cctbx_args}
fi

END_XTC=$(date +"%s")
ELAPSED=$((END_XTC-START_XTC))
echo IOTA_XTC_SingleRank_TimeElapsed ${ELAPSED} ${START_XTC} ${END_XTC} ${RUN}
