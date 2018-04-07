# Building STRUMpack with Spack
Building the dependencies for STRUMpack can be a challenge to perform in a reproducible manner, across many different systems. Using the [Spack](http://spack.readthedocs.io/en/latest/index.html) installer this can easily be used to generate `libstrumpack.a`, and all required dependencies.

To use Spack, the instructions below (be aware, this may take a long time):

```
git clone https://github.com/spack/spack.git && cd spack
./bin/spack create -t generic -f https://github.com/pghysels/STRUMPACK/archive/v2.2.0.tar.gz
```

Spack should now have created a directory structure to hold STRUMpack. Next, copy the file `package.py` into the STRUMpack package directory.

```
cd <SPACK_ROOT>/var/spack/repos/builtin/packages/strumpack/
wget https://raw.githubusercontent.com/ExaFEL/exafel_project/master/95-strumpack_cctbx/spack_installation/package.py
```

This installer is a bare-bones (OpenBLAS, Metis, ParMETIS, SCOTCH, PTSCOTCH, OpenMPI, ScaLAPACK) installation, using the default system compiler. The installation will build all of the above packages, their entire dependency trees, and create a `libstrumpack.a` static library (and headers) in a directory accessible with `./bin/spack location --install-dir strumpack`. All dependencies can be found in a similar manner.

Different MPI versions and compiler versions can be specified, but this requires spending more time with the [Spack documentation](http://spack.readthedocs.io/en/latest/tutorial_packaging.html#creating-the-package-file). It will automatically build any element of the dependency tree, allowing A/B performance testing between compiler versions.
