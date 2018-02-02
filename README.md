# srfpython
  
- programs for surface wave dispersion curves in python
- compute, display, invert 1D depth models
- based on Herrmann codes Computer Program in seismology

install

install tetedenoeud (required)
> move to the installation path (e.g. "~/git") and get tetedenoeud
>
> ```
> cd ~/git
> git clone http://github.com/obsmax/tetedenoeud.git
> ```
>
> follow instructions in tetedenoeud/README.md

install srfpython

> move to the installation path (e.g. "~/git") and get srfpython
>
> ```
> cd ~/git
> git clone http://github.com/obsmax/srfpython.git
> ```
>
> activate the virtual environment created for tetedenoeud
>
> ```
> source activate tetedenoeud
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

add the bin directory to the path (recommended)

> custom the following line and add it to
> your .bashrc or .bash_path (linux) or .profile (mac)
>
> ```
> export PATH=$PATH:"~/git/srfpython/srfpython/bin"
> ```

try the tutorial

> ```
> cd ~/git/srfpython/tutorials/tutorial0
> jupyter notebook 
> ```
