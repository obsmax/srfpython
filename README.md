# srfpython
  
- programs for surface wave dispersion curves in python

install

install tetedoeuf (required)
> move to the installation path : e.g. `cd ~/git`
>
> `git clone http://github.com/obsmax/tetedoeuf.git`
>
> follow instructions in tetedoeuf/README.md

install srfpython

> move to the installation path : e.g. `cd ~/git` 
>  
> `git clone http://github.com/obsmax/srfpython.git`
>
> activate the virtual environment created for tetedoeuf
>
> `source activate tetedoeuf`   
> 
> enter the srfpython directory : e.g. `cd ~/git/srfpython`
>   
> `conda install --yes --file requirements.txt`  
> `pip install -e .`
>
> compile fortran codes
>
> move to the source path for Herrmann codes
>
> e.g. `cd ~/git/srfpython/srfpython/Herrmann/src90`
>
> run the compilation scripts
>
> `bash clean.sh && bash compile.sh`
>
> test
>
> `python ~/git/srfpython/srfpython/Herrmann/Herrmann.py`

add the bin directory to the path

> custom the following line and add it to your bashrc or bash_path   
> `export PATH=$PATH:"~/git/srfpython/srfpython/bin`


