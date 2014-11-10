cp ../../src/cheb/cheb1D.f ../../temp/cheb1D.f
cp ../../src/cheb/evaKer.f ../../temp/evaKer.f
cp ../../src/cheb/checkErr.f ../../temp/checkErr.f

cp ../prini.f ../../temp/prini.f
cp getRankTest.f ../../temp/getRankTest.f
cp ../matcomp.f ../../temp/matcomp.f

cd ../../temp

gfortran getRankTest.f checkErr.f cheb1d.f evaKer.f prini.f\
  matcomp.f -llapack -lblas

./a.out

rm a.out cheb1dtest.f cheb1d.f evaKer.f prini.f matcomp.f
