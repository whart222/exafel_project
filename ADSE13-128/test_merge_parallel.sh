#!/bin/bash

#mpiexec -n 60 libtbx.python /net/dials/raid1/robertb/projects/project1/cctbx/cctbx.xfel/modules/cctbx_project/xfel/merging/application/phil/phil_mpi.py

#libtbx.python /net/dials/raid1/robertb/projects/project1/cctbx/cctbx.xfel/modules/cctbx_project/xfel/merging/application/phil/phil_mpi.py input.path=/net/dials/#raid1/robertb/projects/project1/ld91_int_results/r0095/033/out

#mpiexec -n 60 libtbx.python /net/dials/raid1/robertb/projects/project1/cctbx/cctbx.xfel/modules/cctbx_project/xfel/merging/application/input/file_loader_mpi.py input.path=/net/dials/raid1/robertb/projects/project1/ld91_int_results/r*/033/out

mpiexec -n 60 libtbx.python /net/dials/raid1/robertb/projects/project1/cctbx/cctbx.xfel/modules/cctbx_project/xfel/merging/application/merging/merger_mpi.py input.path=/net/dials/raid1/robertb/projects/project1/ld91_int_results/r0095/033/out

