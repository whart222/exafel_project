#!/bin/bash

#Create most common values for tests
export START=0
export END=1
export EMAIL="loriordan@lbl.gov"

#Set all initial conditions for tests to be the same
for ii in $(seq ${START} ${END}); do
  SDFAC_AUTO[${ii}]="False"
  SDFAC_REFINE[${ii}]="False"
  POSTREF[${ii}]="True"
  ERR_FROM_RES[${ii}]="True"
  AVG_UC[${ii}]="True"
  INCL_NEG[${ii}]="True"
  UC_TOL[${ii}]="0.0065"
  MPI_RANK[${ii}]="50"
  PROG[${ii}]=dev.mpi.cluster_two_merge
done

#Modify values for specific tests
PROG[1]=dev.cxi.mpi_merge_refltable

#Run tests
for ii in 0 ; do
  echo "Started processing ${FILEPREFIX} @ $(date)" | mail -s "$(hostname) PROCESS ${ERR_FROM_RES[$ii]}_${AVG_UC[$ii]}_${INCL_NEG[$ii]}_${UC_TOL[$ii]}_${POSTREF[$ii]}_${SDFAC_AUTO[$ii]}_${SDFAC_REFINE[$ii]}_${MPI_RANK[$ii]}_${PROG[$ii]} STARTED" ${EMAIL}

  export FILEPREFIX=tol${UC_TOL[$ii]}_ErrRes${ERR_FROM_RES[$ii]}_avgUC${AVG_UC[$ii]}_inclNeg${INCL_NEG[$ii]}_postref${POSTREF[$ii]}_sdfacauto${SDFAC_AUTO[$ii]}_sdfacrefine${SDFAC_REFINE[$ii]}_mpirank${MPI_RANK[$ii]}_prog${PROG[$ii]}

  ./do_specific_merge.sh ${FILEPREFIX} ${ERR_FROM_RES[$ii]} ${AVG_UC[$ii]} ${INCL_NEG[$ii]} ${UC_TOL[$ii]} ${POSTREF[$ii]} ${SDFAC_AUTO[$ii]} ${SDFAC_REFINE[$ii]} ${MPI_RANK[$ii]} ${PROG[$ii]} | tee ${FILEPREFIX}.log

  #Create the timing data files for profiling analysis
  cat ${FILEPREFIX}.log | grep -E 'START|END' > ${FILEPREFIX}.dat #Extract timing data
  cut -d ' '  -f1 ${FILEPREFIX}.dat | sort | uniq | sed '/phil/d' > ${FILEPREFIX}_unique.dat #Extract list of unique timing sections, removing a stray #phil
  for u in $(cat ${FILEPREFIX}_unique.dat); do
    echo '' > ${FILEPREFIX}_${u}.dat #empty the file
    cat ${FILEPREFIX}.dat | grep ${u} | cut -f2,4 -d ' ' | sed 's/ TIME=/,/; s/START/0/; s/END/1/' > ${FILEPREFIX}_${u}.dat  #START=0, END=1
  done

  #Format the timing results into a plottable dataset
  for uniq in $(cat ${FILEPREFIX}_unique.dat);do
    echo $uniq
  done

  echo "Finished processing ${FILEPREFIX} @ $(date)" | mail -s "$(hostname) PROCESS ${ERR_FROM_RES[$ii]}_${AVG_UC[$ii]}_${INCL_NEG[$ii]}_${UC_TOL[$ii]}_${POSTREF[$ii]}_${SDFAC_AUTO[$ii]}_${SDFAC_REFINE[$ii]}_${MPI_RANK[$ii]}_${PROG[$ii]} FINISHED" ${EMAIL}

done
