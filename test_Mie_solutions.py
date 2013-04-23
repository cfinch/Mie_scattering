from numpy import *
from bhmie_fortran import bhmie_fortran
from Mie_MATLAB.Mie_MATLAB_interface import call_Mie_MATLAB
from bhmie_matlab.bhmie_matlab_interface import call_bhmie_matlab
import bhmie_herbert_kaiser_july2012

x = 0.099
refrel = 0.75
angles = arange(0.0, pi + 0.01, pi/6)
nang = 7

# ----- Python -----
S1, S2, Qext, Qsca, Qback, gsca = bhmie_herbert_kaiser_july2012.bhmie(x, refrel, angles)
print("bhmie, Python implementation:")
for i in range(len(S1)):
    print("{}   {}  {}".format(cos(angles[i]), S1[i], S2[i]))

# ----- FORTRAN -----
S1, S2, Qext, Qsca, Qback, Gsca = bhmie_fortran.bhmie(x, refrel, nang)
print("bhmie, FORTRAN implementation:")
for i in range(nang):
    print("{}   {}  {}".format(cos(angles[i]), S1[i], S2[i]))

print(Qext)
print(Qback)
print(Gsca)

# ----- MATLAB, one file -----
print('Calling MATLAB single-file implementation')
for a in angles:
    S1, S2 = call_bhmie_matlab(x, refrel, array(a))
    print("{}   {}  {}".format(cos(a), S1, S2))


# ----- MATLAB, multi-file -----
print('Calling MATLAB multi-file')
for a in angles:
    S1, S2 = call_Mie_MATLAB(refrel, x, array(a))
    print("{}   {}  {}".format(cos(a), S1, S2))

