from distutils.core import setup

setup(
    name='srfpython',
    version='0.0',
    packages=['srfpython', 'srfpython.Herrmann', 'srfpython.depthdisp',
              'srfpython.HerrMet', 'srfpython.inversion', 'srfpython.sensitivitykernels'],
    url='',
    license='',
    author='Maximilien Lehujeur',
    author_email='maximilien.lehujeur@gmail.com',
    description='',
    scripts=['srfpython/bin/m96',
             'srfpython/bin/s96',
             'srfpython/bin/HerrMet',
             'srfpython/bin/sker17'])

