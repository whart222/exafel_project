# Notes on using Strumpack within CCTBX
Here are the documented results of using Strumpack on a single node for a variety of data set sizes (`StrumpackSolverMPI_1K`,`StrumpackSolverMPI_5K`,`StrumpackSolverMPI_10K`). All tests were performed on `dials.lbl.gov`, and allow the tests to be repeated at the user's discretion. Example matrices for a variety of different refinement parameters are listed in the given paths, and the times represent a single solution.

# Setting up and running STRUMPACK
To build STRUMPACK alongside a conda cctbx.xfel build, follow the instructions given [here](https://exafel.github.io/docs/psana-cctbx-install)/[here](https://github.com/ExaFEL/exafel_project/tree/master/nks) with the conda packages ammended to the following before building cctbx:

```bash
source activate myEnv;
conda install -y IPython h5py wxpython pillow libtiff mysql-python jinja2 matplotlib scipy mpi4py;
```

Building is now carried out as normal:
```bash
wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py;
python bootstrap.py hot update --builder=xfel --sfuser=<USERNAME> --cciuser=<USERNAME>;
python bootstrap.py build --builder=xfel --with-python=`which python` --nproc=<NUM_CORES>;
source ./build/setpaths.sh;
```

The MPI-enabled work requires the *cctbx_project* git repository be checked-out into the *strumpack_solver_backend* branch (eventually all changes will be merged with upstream *master*).
```bash 
cd modules/cctbx_project
git checkout strumpack_solver_backend
cd ../..
```

We can now proceed with building STRUMPACK, and it's required dependencies. The script located at `github.com/exafel/exafel_project/master/strumpack/STRUMPACK_installer_shared.sh` will acquire all dependencies and build the libraries against the conda MPI installation, or build OpenMPI if this is not available (such as may be true under MacOS).

```bash
wget https://raw.githubusercontent.com/ExaFEL/exafel_project/master/strumpack/STRUMPACK_installer_shared.sh
chmod +x STRUMPACK_installer_shared.sh
./STRUMPACK_installer_shared.sh
```

The STRUMPACK (and dependencies) binmaries, libraries and headers will be installed into `strumpack_build/{bin,lib,include}`, of which will be added to the dispatcher environment given successful completion of the installation script. With the presence of the new libraries, refreshing the dispacther and rebuilding the packaes will alow any STRUMPACK-enabled Boost.Python extension modules to be built.

```bash
cd build && libtbx.refresh
make
```
To test the new libraries several test scripts are available, for a provided A matrix and b vector solution. The solver will compare both the original Eigen-based solver and the new STRUMPACK solvers for both the OpenMP backend (`cctbx_project/scitbx/examples/bevington/strumpack_eigen_solver.py`) and the distributed MPI-enabled solver (`cctbx_project/scitbx/examples/bevington/strumpack_eigen_solver_mpi_dist.py`).

Sample run commands are given in the provided notebooks.

Full integration with Samosa will be provided shortly to allow the solvers to be used therein. Initial work to provide selective choice of the Eigen solver backend is available in the `eigen_solver_algo` branch of `https://github.com/cctbx/cctbx_project/`.

# Scalability tests on Cori

A full-scale test of the scalability of the solvers on Cori for KNL and Haswell with OpenMP for intra-node tests, and MPI+OpeMP for inter-node tests is currently in progress.

# Experimental spack build
An experimental set of commands to build STRUMPACK using spack is given [here](https://github.com/ExaFEL/exafel_project/tree/master/95-strumpack_cctbx/spack_installation). This is not supported, and not guaranteed to work. It was used as an example environment to build STRUMPACK dependencies.
