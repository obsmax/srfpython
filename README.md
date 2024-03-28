# SrfPython
  
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
conda create -n py3-srfpython
conda activate py3-srfpython
python --version  # must be 3.X, or conda failed, retry in new tty
```

- Go to the cloned repository and install the package

```bash
# install in editable mode 
# i.e. changes in python programs do not require re-installing the package
cd ~/git/srfpython
python -m pip install -e .

# test the compilation using 
# (must work anywhere, no need to custom path or anything)
srfpython_tester.py
```

- If you wish to use the jupyter notebooks with python3 (optional):


```bash
# make sure the environment is activated
conda activate py3-srfpython

# install with
conda install decorator=4.4.1 jupyter --yes
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

## Contributing authors
* 2024 : The program has been converted from python 2.7 to >3.7 by Pierric Mora.


## Related references
* Lehujeur, M., Chevrot, S., Villaseñor, A., Masini, E., Saspiturry, N., Lescoutre, R., Sylvander, M., 2021. Three-dimensional shear velocity structure of the Mauléon and Arzacq Basins (Western Pyrenees). BSGF - Earth Sci. Bull. 192, 47. https://doi.org/10.1051/bsgf/2021039
* Lehujeur, M., Vergne, J., Schmittbuhl, J., Zigone, D., Le Chenadec, A., & EstOF Team (2018). Reservoir imaging using ambient noise correlation from a dense seismic network. Journal of Geophysical Research: Solid Earth, 123. https://doi.org/10.1029/2018JB015440

