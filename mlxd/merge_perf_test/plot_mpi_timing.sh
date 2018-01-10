#!/bin/bash
for ii in CORR_EXTEND DMIN_EXTEND REJFRAC_EXTEND WAVELENGTH_EXTEND UCVAL_ADDCELLS DICT_MERGE ISIGI_CID ISIGI_EXTEND CTABLE_EXTEND SEQ_ADD OBS_EXTEND; do
  #Extract the start timing data for the above strings
  cat $1 \
    | grep $ii \
    | tr ';' '\n' \
    | grep START \
    | awk '{ print $3,$4}' \
    | sed 's/RANK=//; s/:/,/; s/ TIME=/,/' \
    > ./REFL_REDUCE_CS_${ii}_START.csv
  #Extract the end timing data for the above strings
  cat $1 \
    | grep $ii \
    | tr ';' '\n' \
    | grep END \
    | awk '{ print $3,$4}' \
    | sed 's/RANK=//; s/:/,/; s/ TIME=/,/' \
    > ./REFL_REDUCE_CS_${ii}_END.csv
  #Pass the extracted data into the python script to generate the timing plot
  libtbx.python plot_mpi_timing.py REFL_REDUCE_CS_${ii}_START.csv REFL_REDUCE_CS_${ii}_END.csv
done
