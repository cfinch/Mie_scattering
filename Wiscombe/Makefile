# Makefile for Warren J. Wiscombe's Mie scattering code
# ftp://climate1.gsfc.nasa.gov/wiscombe/Single_Scatt/Homogen_Sphere/Exact_Mie/

objects = ErrPack.o RDI1MACH.o MIEV0.o MVTstNew.o
target = wiscombe_mie

wiscombe_mie : $(objects)
	gfortran $(objects) -o $@

$(objects) : %.o: %.f
	gfortran -c $< -o $@

clean:
	rm -f wiscombe_mie *.o *.mod
