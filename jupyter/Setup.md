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
Take note of the generated token, and copy it.

Next, set up port forwarding of the Jupyter notebook remote instance to your local machine (run locally):
```bash
ssh -N -f -L localhost:8888:localhost:8889 <remote_username>@<remote_server>
```

Visit `localhost:8888` in your browser, and enter the token from earlier. This should now allow you to run 
your remote Jupyter notebook instance in your local browser.
