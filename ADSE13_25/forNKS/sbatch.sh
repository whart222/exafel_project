# 1. Source DIALS available to you on your computer
DIALS=/net/viper/raid1/asmit/ExaFEL/cctbx.xfel.viper_iota
source ${DIALS}/build/setpaths.sh

# 2. Set IMAGE_PATH i.e where image 64 is located on your computer
IMAGE_PATH=/net/cci/asmit/raid1/iota/LS49_sim/01_stills_process/images_100

# 3. Now run stills_process_modified.py
libtbx.python stills_process_modified.py base_stills.phil ${IMAGE_PATH}/step5_MPIbatch_000064.img.gz
