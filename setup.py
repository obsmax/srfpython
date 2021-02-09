import os, subprocess
import setuptools
from setuptools.command.build_py import build_py
from setuptools.command.develop import develop


# ================ paths and files
version_file = os.path.join('srfpython', 'version.py')
fortran_src_path = os.path.join('srfpython', 'Herrmann', 'src')
packages = setuptools.find_packages()

# checks
assert os.path.isfile(version_file), "{} not found".format(version_file)
assert os.path.isdir(fortran_src_path), "{} not found".format(fortran_src_path)


# ================ get version number from version file
# fuck distutils boolshit, find version string by myself
if not os.path.isfile(version_file):
    raise IOError(version_file)

with open(version_file, "r") as fid:
    for line in fid:
        if line.strip('\n').strip().startswith('__version__'):
            __version__ = line.strip('\n').split('=')[-1].split()[0].strip().strip('"').strip("'")
            break
    else:
        raise Exception('could not detect __version__ affectation in {version_file}'.format(version_file=version_file))

# ================ load description
with open("README.md", "r") as fh:
    long_description = fh.read()


# ============= Custom build_py
def make():
    # run command "make all" in ./srfpython/Herrmann/src
    proc = subprocess.Popen(
        ['/usr/bin/make', 'all'],
        shell=True, cwd=fortran_src_path)
    proc.wait()


class CustomBuilder(build_py):
    # needed if the package is built, i.e. files will be copied to a site-packages directory
    def run(self):
        # it is important to run the make file before building the python packages
        # so that the shared library files (.so) will be copied to the site-packages directory
        # they must appear in the MANIFEST.in file !
        make()
        build_py.run(self)


class CustomDevelop(develop):
    # needed if the package is installed in editable mode

    def run(self):
        make()
        develop.run(self)


setuptools.setup(
    name='srfpython',
    version=__version__,
    packages=packages,
    url='https://github.com/obsmax/srfpython',
    license='',
    author='Maximilien Lehujeur',
    author_email='maximilien.lehujeur@gmail.com',
    description='compute/inverse surface waves dispersion curves, based on Hermmann codes CPS',
    long_description=long_description,
    install_requires=['numpy', 'scipy', 'matplotlib', 'future'],
    python_requires=">=2.7,<2.8",
    classifiers=[
        "Programming Language :: Python :: 2",
        "Operating System :: Linux"],
    cmdclass={"build_py": CustomBuilder,
              "develop": CustomDevelop},
    scripts=['srfpython/bin/m96',
             'srfpython/bin/s96',
             'srfpython/bin/HerrMet',
             'srfpython/bin/sker17',
             'srfpython/bin/Herrmann.py'])
