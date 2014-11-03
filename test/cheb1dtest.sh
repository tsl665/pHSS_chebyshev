cp ../src/cheb/cheb1D.f ../temp/cheb1D.f
cp ../src/cheb/evaKer.f ../temp/evaKer.f

cp prini.f ../temp/prini.f
cp cheb1dtest.f ../temp/cheb1dtest.f
cp matcomp.f ../temp/matcomp.f

cd ../temp

gfortran cheb1dtest.f cheb1d.f evaKer.f prini.f matcomp.f \
	-llapack -lblas

./a.out

rm a.out cheb1dtest.f cheb1d.f evaKer.f prini.f matcomp.f
