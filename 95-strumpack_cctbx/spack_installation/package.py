##############################################################################
# Copyright (c) 2013-2018, Lawrence Livermore National Security, LLC.
# Produced at the Lawrence Livermore National Laboratory.
#
# This file is part of Spack.
# Created by Todd Gamblin, tgamblin@llnl.gov, All rights reserved.
# LLNL-CODE-647188
#
# For details, see https://github.com/spack/spack
# Please also see the NOTICE and LICENSE files for our notice and the LGPL.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License (as
# published by the Free Software Foundation) version 2.1, February 1999.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and
# conditions of the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
##############################################################################
#
# This is a template package file for Spack.  We've put "FIXME"
# next to all the things you'll want to change. Once you've handled
# them, you can save this file and test your package like this:
#
#     spack install strumpack
#
# You can edit this file again by typing:
#
#     spack edit strumpack
#
# See the Spack documentation for more information on packaging.
# If you submit this package back to Spack as a pull request,
# please first remove this boilerplate and all FIXME comments.
#
from spack import *


class Strumpack(CMakePackage):
    """STRUctured Matrix PACKage"""

    # FIXME: Add a proper url for your package's homepage here.
    homepage = "http://portal.nersc.gov/project/sparse/strumpack/"
    url      = "https://github.com/pghysels/STRUMPACK/archive/v2.2.0.tar.gz"

    version('2.2.0', 'f1174de1c70da0c209ab19e19d80eff3')

    # FIXME: Add dependencies if required.
    depends_on('openmpi',type=('build','run','link'))
    depends_on('scotch', type=('build', 'run','link'))
    depends_on('parmetis', type=('build','run','link'))
    depends_on('openblas', type=('build','run','link'))
    depends_on('scalapack', type=('build','run','link'))
    patch('cmakelist_ctest.patch', level=1, working_dir='.')
    
    def cmake_args(self):
      spec = self.spec
      options = []
      options.extend([
          '-DCMAKE_INSTALL_PREFIX=%s' % (self.prefix), 
          '-DCMAKE_POSITION_INDEPENDENT_CODE:BOOL=true',
          '-DBUILD_TESTING=%s' % ('OFF'), 
          '-DSTRUMPACK_DEV_TESTING=%s' % ('OFF'), 
          '-DCMAKE_BUILD_TYPE=%s' % ('Debug'), 
          '-DSTRUMPACK_USE_OPENMP=ON', 
          '-DSTRUMPACK_USE_PARMETIS=ON', 
          '-DSTRUMPACK_USE_SCOTCH=ON',
          '-DCMAKE_CXX_FLAGS=%s %s %s %s' % ('-Wall', '-Wfatal-errors', '-Wextra', '-Wno-unused-parameter'), 
          '-DMETIS_INCLUDES=%s' % spec['metis'].prefix.include,
          '-DMETIS_LIBRARIES=%s' % join_path( spec['metis'].prefix.lib, 'libmetis.so'), 
          '-DPARMETIS_INCLUDES=%s' % spec['parmetis'].prefix.include,
          '-DPARMETIS_LIBRARIES=%s' % join_path(spec['parmetis'].prefix.lib, 'libparmetis.so'), 
          '-DSCALAPACK_LIBRARIES=%s' % join_path( spec['scalapack'].prefix.lib, 'libscalapack.so' ),
          '-DLAPACK_LIBRARIES=%s' % join_path( spec['openblas'].prefix.lib, 'libopenblas.so' ), 
          '-DBLAS_LIBRARIES=%s' % join_path( spec['openblas'].prefix.lib, 'libopenblas.so' ), 
          '-DSCOTCH_INCLUDES=%s' % spec['scotch'].prefix.include, 
          '-DSCOTCH_LIBRARIES="%s;%s;%s;%s"' % 
          ( join_path(spec['scotch'].prefix.lib, 'libscotch.so'), 
            join_path(spec['scotch'].prefix.lib, 'libscotcherr.so'),
            join_path(spec['scotch'].prefix.lib, 'libptscotch.so'),
            join_path(spec['scotch'].prefix.lib, 'libptscotcherr.so') 
          )
        ])
      return options

    def setup_environment(self, spack_env, run_env):
      spec = self.spec

      spack_env.set('CC', spec['mpi'].mpicc)
      spack_env.set('FC', spec['mpi'].mpifc)
      spack_env.set('CXX', spec['mpi'].mpicxx)

#    def install(self, spec, prefix):
#        make()
#        make('install')
