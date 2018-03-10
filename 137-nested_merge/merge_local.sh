#!/bin/bash
export trial=${TAG}
cd ${MERGE_ROOT}/${TAG}

export effective_params="d_min=2.0 \
targlob=\"${TARDATA}\" \
model=${MERGE_ROOT}/4ngz.pdb \
backend=FS \
mysql.runtag=${trial} \
scaling.report_ML=True \
pixel_size=0.11 \
nproc=1 \
postrefinement.enable=True \
scaling.mtz_file=${MERGE_ROOT}/4ngz.mtz \
scaling.mtz_column_F=f(+) \
min_corr=-1.0 \
unit_cell_length_tolerance=0.0065 \
raw_data.errors_from_sample_residuals=False \
raw_data.sdfac_refine=True \
raw_data.error_models.sdfac_refine.random_seed=42 \
cell_rejection.unit_cell=79.0,79.0,38.1,90,90,90 \
cell_rejection.space_group=P43212 \
output.prefix=${trial}"

#source ${CONDA_ROOT}/build/setpaths.sh
if [ ${MULTINODE} == "True" ]; then
  # source the cctbx build
  echo ${effective_params}
  dev.mpi.cluster_two_merge ${effective_params}

  #To manually specify the merge file, use the following format
  #libtbx.python my_merge.py ${effective_params}
  exit
fi

#source the phenix build for single node execution:
cxi.xmerge ${effective_params} && \

source ${PHENIX_ROOT}/phenix-installer-dev-2880-source/build/setpaths.sh 

phenix.xtriage ${trial}.mtz scaling.input.xray_data.obs_labels=Iobs |tee xtriage_${trial}.out && \

phenix.refine ${EXAFEL_DIR}/exafel_project/nks/Refine_4/LM14_1colorYblyso_refine_4.pdb \
refinement.output.prefix=${trial} \
xray_data.file_name=${trial}.mtz \
xray_data.labels="Iobs" \
xray_data.r_free_flags.file_name=${EXAFEL_DIR}/exafel_project/nks/merge118/reflections.mtz \
xray_data.r_free_flags.label=FreeR_flag \
main.number_of_macro_cycles=3 optimize_xyz_weight=True \
optimize_adp_weight=True nproc=20 \
ordered_solvent=True ordered_solvent.mode=every_macro_cycle \
refine.adp.individual.anisotropic="element Yb" \
refinement.input.symmetry_safety_check=warning \
refine.strategy="*individual_sites *individual_sites_real_space *rigid_body *individual_adp group_adp tls *occupancies group_anomalous" \
xray_data.force_anomalous_flag_to_be_equal_to=True \
main.wavelength=1.3853 --overwrite && \

libtbx.python ${EXAFEL_DIR}/exafel_project/nks/map_height_at_atoms.py \
${trial}_001.pdb \
${trial}.mtz \
input.xray_data.labels=Iobs \
input.symmetry_safety_check=warning \
xray_data.r_free_flags.file_name=${trial}_001.mtz \
xray_data.r_free_flags.label=R-free-flags \
selection="element Yb" |tee ${trial}_peakht.log && \

phenix.molprobity input.pdb.file_name=${trial}_001.pdb \
molprobity.flags.clashscore=False molprobity.flags.ramalyze=False \
molprobity.flags.omegalyze=False molprobity.flags.rotalyze=False molprobity.flags.cbetadev=False \
molprobity.flags.nqh=False molprobity.flags.rna=False molprobity.flags.restraints=False \
output.coot=False output.probe_dots=False output.prefix=${trial}_molprobity && \

libtbx.python ${EXAFEL_DIR}/exafel_project/nks/json/to_json.py ${MERGE_ROOT} ${trial}

