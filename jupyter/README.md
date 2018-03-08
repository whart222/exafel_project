# Remote server invocation
To use Jupyter remotely with cctbx, first install Jupyter notebook 
on the remote machine where the computation will run (assuming conda installation):
```bash
conda install jupyter
```

Next, launch a Jupyter notebook instance using the `libtbx.ipython` dispatcher to allow all access to all
cctbx functionality:
```bash
libtbx.ipython notebook --no-browser --port=8889
```
The above command assumes that ipython was installed with the conda environment. If it was not, install it with `conda instal ipython` and then run `libtbx.refresh` to update the dispatcher to generate `libtbx.ipython`.
Take note of the generated token by running the notebook command, and copy it.

Next, set up port forwarding of the Jupyter notebook remote instance to your local machine (run locally):
```bash
ssh -N -f -L localhost:8888:localhost:8889 <remote_username>@<remote_server>
```
Visit `localhost:8888` in your browser, and enter the token from earlier. This should now allow you to run 
your remote Jupyter notebook instance in your local browser.

# NERSC Cori invocation
To use Jupyter notebooks on Cori, one can visit https://jupyter-dev.nersc.gov where logging in presents the Jupyter dashboard on the user's home directory. From here, it is possible to start a notebook with a custom Python distribution, and still have access to the cluster queing system.

Assuming that we have a conda build of cctbx, we can configure Jupyter to use *libtbx.python* as a kernel for executing the notebook. We first create a new directory to store out kernel information with 
```bash
mkdir $HOME/.local/share/jupyter/kernels/kernel_name/
```
Next, we create a JSON file configuring the required kernel paths, named `kernel.json`, and place it into the above directory. The Structure of the file should be as follows:
```json
{
 "display_name": "libtbx.python",
 "language": "python",
 "argv": [
"/global/cscratch1/sd/mlxd/libtbx_env.sh",
  "-f",
  "{connection_file}"
 ]
}
```
where I have created a helper script, `libtbx_env.sh`, to set the environment as:
```bash
#!/usr/bin/env bash
source <PATH_TO_CONDA_DIR>/bin/activate myEnv
source <PATH_TO_CCTBX_DIR>/build/setpaths.sh

exec libtbx.python -m ipykernel $@
```

With these files in place, the Jupyter hub dashboard should allow the new kernel `libtbx.python` to be started on a Cori server.
It can be useful also to set up integration with the SLURM environment within the notebook. We can use the NERSC Github repo https://github.com/NERSC/slurm-magic to allow for this. 

From a Cori shell, with libtbx.python available, we first build and install the package into the Python distribution.
```
git clone https://github.com/NERSC/slurm-magic
libtbx.python -m pip install ./slurm-magic
conda install pandas -y
```
The above packagaes allow us to directly access and submit jobs to the NERSC SLURM queueing system. These commands take the form of IPython *magic* commands, and are accessible as `%command`. Here, we must load the new magic package, then register it with the IPython environment.
```python
from slurm_magic import SlurmMagics
ip = get_ipython()
ip.register_magics(SlurmMagics)
```
We may now access the full breadth of SLURM commands, such as `%sbatch` to submit jobs, `%squeue` to examine the queue. A full list of the supported commands are available on the slurm-magic Github repo.
