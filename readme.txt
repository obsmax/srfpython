srfpython
Maximilien Lehujeur, 18/04/2017
maximilien.lehujeur@gmail.com

Compute surface wave dispersion curves for Love and Rayleigh waves in python2.7
the heart of the code is in fortran77 modified after Computer Program in Seismology, Herrmann and Ammon, 2002
the original codes (srfpre96, srfdis96) were slightly modified to prevent disk ios and communicate through stdin and stdout

dependencies 
    numpy 
    (matplotib)

I/ install (linux)

1) compile fortran codes, needs gcc gfortran
>> cd src
>> bash clean.sh
>> bash compile.sh
>> ls ../bin   

1.1) optional : test the fortran codes with srfpyhon/src/test/test.sh

2) add srfpython/bin to your python path

3) in python, import function srfdis17 from module srfid17
>> from srfdis17 import srfdis17
>> help(srfdis17)

II/ run demo (see __main__ section in srfpython/bin/srfdis17)

cd srfpython/bin
python srfdis17.py
