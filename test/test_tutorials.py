"""
Script used to convert all files from python2 to python 3, Pierric Mora 2024
Walks the srfpython/ directory and executes all .py files
"""

import os
import subprocess

tuto_dir = '../tutorials'

to_test = ['00_simple_dispersion_example/00_using_scripts/main_script.sh',
           '00_simple_dispersion_example/01_using_python_programs/000_create_dephmodel.py',
           '00_simple_dispersion_example/01_using_python_programs/001_forward_dispersion.py',
           '00_simple_dispersion_example/02_1d_sensitivity_kernel.py',
           '01_simple_inversion_example/main_script.sh',
           '03_cube_inversion_example/main_script.sh',
           '10_notebooks/00_how_to_use_srfpython.ipynb'
           ]

ignore = ('03_cube_inversion_example', )

to_test = filter(lambda s: all((not(skip in s) for skip in ignore)), to_test)

CWD = os.getcwd()

for f in to_test:
    fpath = os.path.join(tuto_dir, f)
    fdir = os.path.dirname(fpath)
    fname = os.path.basename(fpath)
    os.chdir(fdir)
    print(f'##########\n#####\n# executing {fpath}...')
    if fpath.endswith('.py'):
        r = subprocess.run(('python3', fname))
    elif fpath.endswith('.sh'):
        r = subprocess.run(('bash', fname))
    elif fpath.endswith('.ipynb'):
        r = subprocess.run(('jupyter', 'nbconvert', '--execute', fname, '--to', 'html'))
    r.check_returncode()
    os.chdir(CWD)

# clear outpus
# subprocess.run(('rm', '-f', '000.mod96', '_HerrMet.param', 'toto.bdd'))
