# Docker/Shifter
The included Dockerfile will generate a container of the cctbx.xfel build. 
To build, follow the installation directions to acquire the Miniconda installer and the cctbx modules directory. Next, the MPICH sources and MPI4PY sources should also be acquired. 
```bash
curl -O http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz
curl -L -O https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-3.0.0.tar.gz
``` 

Running `docker build .` will then build the container with the latest available sources in the current directory.

To upload the image to NERSC, please following the directions listed [here](http://www.nersc.gov/users/software/using-shifter-and-docker/using-shifter-at-nersc/).
