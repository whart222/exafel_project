#!/bin/bash
START_MERGE=$(date +"%s")
TRIAL=0

# LD91 data directory
LD91_DATA_DIR=

MERGE_ROOT=$PWD

# output directory
OUT_DIR=$PWD

# trial
TRIAL_F="$(printf "%03d" ${TRIAL})"

# setup playground
mkdir -p ${OUT_DIR}/${TRIAL_F}/out
mkdir -p ${OUT_DIR}/${TRIAL_F}/stdout
mkdir -p ${OUT_DIR}/${TRIAL_F}/tmp

export effective_params_nks_ex="\
input.path=${LD91_DATA_DIR}/r*/033/out \
input.experiments_suffix=_integrated_experiments.json \
input.reflections_suffix=_integrated.pickle \
input.parallel_file_load.method=uniform \
filter.algorithm=unit_cell \
filter.unit_cell.value.target_unit_cell=79.0,79.0,38.3,90,90,90 \
filter.unit_cell.value.target_space_group=P43212 \
filter.unit_cell.value.relative_length_tolerance=0.02 \
filter.outlier.min_corr=-1.0 \
select.algorithm=significance_filter \
scaling.unit_cell=78.589,78.589,37.017,90,90,90 \
scaling.space_group=P43212 \
scaling.model=${MERGE_ROOT}/4ngz.pdb \
scaling.mtz.mtz_column_F=f(+) \
scaling.resolution_scalar=0.96 \
postrefinement.enable=True \
postrefinement.algorithm=rs \
merging.d_min=2.0 \
merging.merge_anomalous=False \
merging.set_average_unit_cell=True \
merging.error.model=errors_from_sample_residuals \
statistics.n_bins=10 \
statistics.cciso.mtz_file=${MERGE_ROOT}/4ngz.mtz \
statistics.cciso.mtz_column_F=f(+) \
output.prefix=${TRIAL_F} \
output.output_dir=${OUT_DIR}/${TRIAL_F}/out \
output.tmp_dir=${OUT_DIR}/${TRIAL_F}/tmp \
output.do_timing=True \
output.log_level=1 \
parallel.a2a=1"

#run merge
cctbx.xfel.merge ${effective_params_nks_ex}

END_MERGE=$(date +"%s")
ELAPSED=$((END_MERGE-START_MERGE))
echo TotalElapsed_OneCore ${ELAPSED} ${START_MERGE} ${END_MERGE}

