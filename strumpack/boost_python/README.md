# Build instructions on Cori

- Ensure the GCC programming environment is loaded, and gcc4.9.3 is selected.
```bash
module swap PrgEnv-intel PrgEnv-gnu; 
module swap gcc gcc/4.9.3
```
- Add the built STRUMPACK library path to LD_LIBRARY_PATH
```bash
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/global/cscratch1/sd/mlxd/feb_sprint/modules/strumpack/builds/lib/
```
- Compile the Boost.Python extension module against a pre-existing xfel build of cctbx
```bash
g++ -c -fPIC ./strumpack_solver_ext.cc -I/global/cscratch1/sd/mlxd/feb_sprint/modules/strumpack/builds/include/ -I$(cc --cray-print-opts=cflags) -std=c++11 -D_GLIBCXX_USE_CXX11_ABI=0 -DBOOST_PYTHON_MAX_BASES=2 -I/global/cscratch1/sd/mlxd/feb_sprint/modules/boost -I/global/cscratch1/sd/mlxd/feb_sprint/modules/boost -fPIC -fno-strict-aliasing -D_GLIBCXX_USE_CXX11_ABI=0 -Wall -Wno-sign-compare -Wno-unknown-pragmas -Wno-parentheses -Winit-self -Wno-unused-local-typedefs -Werror=vla -DNDEBUG -O3 -funroll-loops -ffast-math -DBOOST_ALL_NO_LIB -I/global/cscratch1/sd/mlxd/feb_sprint/build/include -I/global/cscratch1/sd/mlxd/feb_sprint/modules/cctbx_project -I/global/cscratch1/sd/mlxd/merge_perf/miniconda/envs/strumpack/include/python2.7 -I/global/cscratch1/sd/mlxd/feb_sprint/modules -I/global/cscratch1/sd/mlxd/feb_sprint/modules/cctbx_project -o strumpack_solver_ext.o
```
- Link the libraries and produce `strumpack_solver.so`, which can now be loaded into python:
```bash
g++ -o strumpack_solver.so -shared -s strumpack_solver_ext.o -Llib -L/global/cscratch1/sd/mlxd/feb_sprint/modules/strumpack/builds/lib/  -L/global/cscratch1/sd/mlxd/feb_sprint/build/lib -L/global/cscratch1/sd/mlxd/feb_sprint/modules/cctbx_project/lib -lboost_python -lm -lscitbx_boost_python -lboost_python -lcctbx -lstrumpack
```

```python
from sparse_solver_ext import sparse_solver as ss
res = ss(5)
# Initializing STRUMPACK
# running serially, no OpenMP support!
# number of tasking levels = 9 = log_2(#threads) + 3
# initial matrix:
#   - number of unknowns = 25
#   - number of nonzeros = 105
# nested dissection reordering:
#   - Geometric reordering
#   - strategy parameter = 8
#   - number of separators = 7
#   - number of levels = 3
#   - nd time = 9.4405e-05
#   - symmetrization time = 3.596e-06
# symbolic factorization:
#   - nr of dense Frontal matrices = 7
#   - nr of HSS Frontal matrices = 0
#   - symb-factor time = 1.2437e-05
# multifrontal factorization:
#   - factor time = 1.05475
#   - factor nonzeros = 265
#   - factor memory = 0.00212 MB
#   - factor memory/nonzeros = 100 % of multifrontal
#   - maximum HSS rank = 0
#   - HSS compression = false
#   - relative compression tolerance = 0.01
#   - absolute compression tolerance = 1e-08
#   - normal(0,1) distribution with minstd_rand engine
REFINEMENT it. 0  res =            5  rel.res =            1  bw.error =            1
REFINEMENT it. 1  res =  3.62144e-15  rel.res =  7.24288e-16  bw.error =  1.80411e-16
# DIRECT/GMRES solve:
#   - abs_tol = 1e-10, rel_tol = 1e-06, restart = 30, maxit = 5000
#   - number of Krylov iterations = 1
#   - solve time = 0.0296811
# COMPONENTWISE SCALED RESIDUAL = 1.6237e-16
In [6]: res.b
TypeError: No to_python (by-value) converter found for C++ type: std::vector<double, std::allocator<double> >
```
- Vectors `b` and `x` can be directly accessed after solution. Modifications to ensure Python can read these objects will be required.
