from distutils.core import setup
from srfpython.version import __version__

setup(
    name='srfpython',
    version=__version__,
    packages=['srfpython', 'srfpython.Herrmann', 'srfpython.depthdisp',
              'srfpython.HerrMet', 'srfpython.inversion', 'srfpython.sensitivitykernels'],
    url='',
    license='',
    author='Maximilien Lehujeur',
    author_email='maximilien.lehujeur@gmail.com',
    description='',
    scripts=["./srfpython/bin/m96",
             "./srfpython/bin/s96",
             "./srfpython/bin/sker17",
             "./srfpython/bin/HerrMet"])
