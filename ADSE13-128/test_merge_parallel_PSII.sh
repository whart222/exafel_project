#!/bin/bash
START_MERGE=$(date +"%s")

TRIAL=0

#cctbx
#source /build/setpaths.sh

# output directory
OUT_DIR=${SCRATCH}/project1/work4

# trial
TRIAL_F="$(printf "%03d" ${TRIAL})"

# setup playground
mkdir -p ${OUT_DIR}/${TRIAL_F}/out
mkdir -p ${OUT_DIR}/${TRIAL_F}/stdout
mkdir -p ${OUT_DIR}/${TRIAL_F}/tmp

libtbx.python /global/cscratch1/sd/robbol/project1/modules/cctbx_project/xfel/merging/application/merging/merger_mpi.py input.path=/global/cscratch1/sd/robbol/PSII/*/r*/combine_experiments_*/intermediates input.reflections_suffix=_reintegrated_reflections.pickle input.experiments_suffix=_reintegrated_experiments.json filter.unit_cell.value.target_unit_cell=117.87,223.14,310.71,90,90,90 filter.unit_cell.value.target_space_group=P212121 filter.resolution.d_min=2.0 output.output_dir=${OUT_DIR}/${TRIAL_F}/out merging.merge_anomalous=False

END_MERGE=$(date +"%s")
ELAPSED=$((END_MERGE-START_MERGE))
echo TotalElapsed_OneCore ${ELAPSED} ${START_MERGE} ${END_MERGE}

