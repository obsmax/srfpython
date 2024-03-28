"""
Script used to convert all files from python2 to python 3, Pierric Mora 2024
Walks the srfpython/ directory and executes all .py files
"""
raise Exception('security (A lot of files need to be cleaned up first), TODO: move to pytest')
import os
import subprocess

srfpython_dir = '../srfpython/'
ignore = ('testbox/tester.py', 'synthetics/synthetics.py')

for root, dirs, files in os.walk(srfpython_dir):
    pyfiles = filter(lambda s: s.endswith('.py'), files)
    pyfiles = filter(lambda s: all((not(s in skipped) for skipped in ignore)), pyfiles)
    for f in pyfiles:
        fpath = os.path.join(root, f)
        print(f'##########\n#####\n# executing {fpath}...')
        r = subprocess.run(('python3', fpath))
        r.check_returncode()

# clear outpus
subprocess.run(('rm', '-f', '000.mod96', '_HerrMet.param', 'toto.bdd'))
