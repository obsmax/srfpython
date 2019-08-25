# srfpython
  
- programs for surface wave dispersion curves in python
- compute, display, invert 1D depth models
- based on Herrmann codes Computer Program in seismology

if you use this program, reference would be greatly appreciated  
how to cite :

> Lehujeur, M., Vergne, J., Schmittbuhl, J., Zigone, D., Le Chenadec, A., & EstOF Team (2018). Reservoir imaging using ambient noise correlation from a dense seismic network. Journal of Geophysical Research: Solid Earth, 123. https://doi.org/10.1029/2018JB015440

install

> move to the installation path (e.g. "~/git") and get srfpython
>
> ```
> cd ~/git
> git clone http://github.com/obsmax/srfpython.git
> ```
>
> create the virtual environment and activate it
>
> ```
> conda create -n srfpython python=2.7
> conda activate srfpython
> # source activate srfpython # on old versions of anaconda
> ```
>
> move to the repository, install the requirements and install the package
>
> ```
> cd ~/git/srfpython
> conda install --yes --file requirements.txt
> pip install -e .
> ```
>
> compile fortran codes
>
> ```
> cd ~/git/srfpython/srfpython/Herrmann/src90
> ./clean.sh 
> ./compile.sh
> ```
>
> test fortran codes using
>
> ```
> python ~/git/srfpython/srfpython/Herrmann/Herrmann.py
> ```

if you plan to use jupyter notebooks with python2 (optional)

> make sure the environment is activated
> ```
> source activate srfpython
> ```
>
> install with
>
> ```
> conda install --yes notebook ipykernel
> ipython kernel install --user
> ```
>

try the tutorials

> ```
> cd ~/git/srfpython/tutorials/tutorial0
> jupyter notebook 
> ```
