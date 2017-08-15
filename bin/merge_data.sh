#!/bin/bash

#This file, while not fully necessary with cxi.merge on the path, is used to showcase the sample parameters for processing cxid9114 with a single KNL node. The model file can be added manually from the PDB, or downloaded using phenix.model 

#DATA_DIR contains the path to the integrated pickle files that will be used in merging
DATA_DIR=/global/cscratch1/sd/mlxd/DataProc/cxid9114/processing/batch_metrology_r0113_013/r0113/000/out/untar/r0113/000/out/

NPROC=68 # Total number of cores (single node) to process using
OUT_FILE=cxid9114_out.mtz #Output file
MODEL=4ngz.pdb #Model used for phase information
PREFIX=cxid9114_merged #String prepended to output files
D_MIN=4.0 #Minimum feature size

#Perform merging step
cxi.merge d_min=${D_MIN} data=${DATA_DIR} pixel_size=0.11 output.prefix=${PREFIX} nproc=${NPROC} backend=FS scaling.mtz_file=${OUT_FILE} -v -e 10 model=${MODEL}
