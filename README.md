# srfpython
  
- programs for surface wave dispersion curves in python
- compute, display, invert 1D depth models
- based on Herrmann codes Computer Program in seismology


## Install

### Pre-requisties
* git
* Anaconda3
* gfortran

### Install instructions
- Go to the installation path (e.g. "~/git") and get srfpython
```bash
cd ~/git
git clone http://github.com/obsmax/srfpython.git
```

- Create the virtual environment and activate it

```bash
conda create -n py27-srfpython python=2.7 --yes
conda activate py27-srfpython
python --version  # must be 2.7.X, or conda failed, retry in new tty
```

- Go to the cloned repository and install the package

```bash
# install in editable mode 
# i.e. changes in python programs do not require re-installing the package
cd ~/git/srfpython
python -m pip install -e .

# test the compilation using 
srfpython_tester.py
```

- If you wish to use the jupyter notebooks with python2 (optional):


```bash
# make sure the environment is activated
conda activate py27-srfpython

# install with
conda install --yes notebook ipykernel
ipython kernel install --user
```

## Tutorials
for a simple dispersion curve example 
```bash
cd ~/git/srfpython/tutorials/00_simple_dispersion_example
# follow instructions in readme.txt 
```

for an depth inversion example 
```bash
cd ~/git/srfpython/tutorials/01_simple_inversion_example
# follow instructions in readme.txt 
```

more in the jupyter notebooks
```bash
cd ~/git/srfpython/tutorials/10_notebooks
jupyter notebook 
```

## Related references

* Lehujeur, M., Vergne, J., Schmittbuhl, J., Zigone, D., Le Chenadec, A., & EstOF Team (2018). Reservoir imaging using ambient noise correlation from a dense seismic network. Journal of Geophysical Research: Solid Earth, 123. https://doi.org/10.1029/2018JB015440
* ...