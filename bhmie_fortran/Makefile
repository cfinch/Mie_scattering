callbhmie : callbhmie.o bhmie.o
	gfortran callbhmie.o bhmie.o -o callbhmie

callbhmie.o : callbhmie.f
	gfortran -c callbhmie.f

bhmie.o : bhmie.f
	gfortran -c bhmie.f

python: bhmie.f
	f2py2.7 -c -m bhmie_fortran bhmie.f 

clean :
	rm -f callbhmie.o bhmie.o callbhmie callbhmie.out *.pyc bhmie_fortran.so
