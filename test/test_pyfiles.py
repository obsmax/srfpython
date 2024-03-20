"""
Walks the srfpython/ directory and executes all .py files
"""

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
