export CCTBX_XFEL=~/cctbx.xfel
export PHENIX_ROOT=~/phenix
# Taken from Exafel/exafel_project/nks
export tag="cxic0415"
export n=$1
export trial=${tag}_${n}
# Set Base Folder where processing and merging is done. Folder structure should be
# BASE_PATH --> merging --> reference
#           --> r0050 ---> XXX_rgYYY --> out
export BASE_PATH=/global/cscratch1/sd/asmit/iota_demo/cxic0415 
###export BASE_PATH=/project/projectdirs/lcls/asmit/iota_demo/cxic0415 
source ${CCTBX_XFEL}/build/setpaths.sh

export effective_params="d_min=1.4 \
output.n_bins=20 \
model=${BASE_PATH}/merging/reference/5jd2.pdb \
targlob=${BASE_PATH}/r0050/014_rg002/out/int-x-${n}.pickle.tar \
backend=FS \
nproc=68 \
scaling.mtz_file=${BASE_PATH}/merging/reference/5jd2-sf.mtz \
scaling.mtz_column_F=i(+) \
postrefinement.enable=True \
postrefinement.algorithm=rs \
min_corr=-1.0 \
include_negatives=True \
set_average_unit_cell=True \
raw_data.sdfac_auto=True \
output.prefix=${trial}"


# Run cxi.merge and cxi.xmerge
cxi.merge ${effective_params}
