#!/bin/bash
START_MERGE=$(date +"%s")
TRIAL=0

# PSII data directory
PSII_DATA_DIR=

MERGE_ROOT=$PWD

# output directory
OUT_DIR=$PWD

# trial
TRIAL_F="$(printf "%03d" ${TRIAL})"

# setup playground
mkdir -p ${OUT_DIR}/${TRIAL_F}/out
mkdir -p ${OUT_DIR}/${TRIAL_F}/stdout
mkdir -p ${OUT_DIR}/${TRIAL_F}/tmp

export input_params="\
input.path=${PSII_DATA_DIR}/m*/r*/combine_experiments_*/intermediates \
input.reflections_suffix=reintegrated_reflections.pickle \
input.experiments_suffix=reintegrated_experiments.json \
input.parallel_file_load.method=node_memory \
input.parallel_file_load.node_memory.architecture=summit \
input.parallel_file_load.node_memory.limit=480.0 \
input.parallel_file_load.node_memory.pickle_to_memory=3.5 \
input.parallel_file_load.ranks_per_node=42 \
input.parallel_file_load.balance=global \
input.parallel_file_load.balance_mpi_alltoall_slices=2 \
filter.algorithm=unit_cell \
filter.unit_cell.value.target_unit_cell=117.0,223.0,307.9,90,90,90 \
filter.unit_cell.value.target_space_group=P212121 \
filter.unit_cell.value.relative_length_tolerance=0.01 \
filter.outlier.min_corr=-1.0 \
select.algorithm=significance_filter \
select.significance_filter.sigma=0.1 \
scaling.unit_cell=117.0,223.0,307.9,90,90,90 \
scaling.space_group=P212121 \
scaling.model=${MERGE_ROOT}/LS11_LS34_LQ39_LN84_LM51_all_OEC_1.92_1022_30.pdb \
scaling.mtz.mtz_column_F=i-obs \
scaling.resolution_scalar=0.96 \
postrefinement.enable=True \
postrefinement.algorithm=rs \
merging.d_min=2.0 \
merging.merge_anomalous=False \
merging.set_average_unit_cell=True \
merging.error.model=ha14 \
statistics.n_bins=20 \
statistics.cciso.mtz_file=${MERGE_ROOT}/noanom_PSII_LM51_LN84_LQ39_LS34_LS10_all_flash_states_0d03_s0_mark0.mtz \
statistics.cciso.mtz_column_F=i-obs \
output.prefix=${TRIAL_F} \
output.output_dir=${OUT_DIR}/${TRIAL_F}/out 
output.tmp_dir=${OUT_DIR}/${TRIAL_F}/tmp \
output.do_timing=True \
output.log_level=1 \
parallel.a2a=1"

#run merge
cctbx.xfel.merge ${input_params}

END_MERGE=$(date +"%s")
ELAPSED=$((END_MERGE-START_MERGE))
echo TotalElapsed_OneCore ${ELAPSED} ${START_MERGE} ${END_MERGE}

