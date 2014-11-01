cp ../src/cheb1D.f ../temp/cheb1D.f
cp ../src/evaKer.f ../temp/evaKer.f

cp prini.f ../temp/prini.f
cp cheb1dtest.f ../temp/cheb1dtest.f

cd ../temp

gfortran cheb1dtest.f cheb1d.f evaKer.f prini.f \
	-llapack -lblas

./a.out

rm a.out cheb1dtest.f cheb1d.f evaKer.f prini.f 
