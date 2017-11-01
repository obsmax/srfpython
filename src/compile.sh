gfortran -O -c max_igetmod.f
gfortran -O -c max_iputmod.f
gfortran -O -c max_lgstr.f
gfortran -O  -c max_setdsp.f
gfortran -O  -c max_setmod.f
gfortran -O  -c max_srfpre96.f
gfortran -O -c max_srfdis96.f


gfortran -O  -o max_srfpre96 max_srfpre96.o max_lgstr.o \
	max_igetmod.o max_iputmod.o max_setdsp.o max_setmod.o 


gfortran -O -o max_srfdis96 max_srfdis96.o  max_igetmod.o 

rm -f max_igetmod.o
rm -f max_iputmod.o
rm -f max_lgstr.o
rm -f max_setdsp.o
rm -f max_setmod.o
rm -f max_srfdis96.o
rm -f max_srfpre96.o

mkdir ../bin
mv max_srfdis96 max_srfpre96 ../bin
cd ../bin
ln -s ../src/srfdis17.py
ln -s ../src/sker17.py
touch __init__.py

ls -l --color ../bin
