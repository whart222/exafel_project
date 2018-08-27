#!/bin/bash

#mpiexec -n 60 libtbx.python /net/dials/raid1/robertb/projects/project1/cctbx/cctbx.xfel/modules/cctbx_project/xfel/merging/application/phil/phil_mpi.py

#libtbx.python /net/dials/raid1/robertb/projects/project1/cctbx/cctbx.xfel/modules/cctbx_project/xfel/merging/application/phil/phil_mpi.py input.path=/net/dials/#raid1/robertb/projects/project1/ld91_int_results/r0095/033/out

#mpiexec -n 60 libtbx.python /net/dials/raid1/robertb/projects/project1/cctbx/cctbx.xfel/modules/cctbx_project/xfel/merging/application/input/file_loader_mpi.py input.path=/net/dials/raid1/robertb/projects/project1/ld91_int_results/r*/033/out

#mpiexec -n 4 libtbx.python /net/dials/raid1/robertb/projects/project1/cctbx/cctbx.xfel/modules/cctbx_project/xfel/merging/application/merging/merger_mpi.py input.path=/net/dials/raid1/robertb/projects/project1/ld91_int_results/r0095/033/out filter.unit_cell.value.target_unit_cell=78.9,78.9,38.1,90,90,90 filter.unit_cell.value.target_space_group=P43212 filter.resolution.d_min=1.7

mpiexec -n 60 libtbx.python /net/dials/raid1/robertb/projects/project1/cctbx/cctbx.xfel/modules/cctbx_project/xfel/merging/application/merging/merger_mpi.py input.path=/net/dials/raid1/robertb/projects/project1/ld91_int_results/r*/033/out filter.unit_cell.value.target_unit_cell=78.9,78.9,38.1,90,90,90 filter.unit_cell.value.target_space_group=P43212 filter.resolution.d_min=1.7

