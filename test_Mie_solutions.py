from numpy import *
from bhmie_fortran import bhmie_fortran
from Mie_MATLAB.Mie_MATLAB_interface import call_Mie_MATLAB
from bhmie_matlab.bhmie_matlab_interface import call_bhmie_matlab
import bhmie_herbert_kaiser_july2012

x = 0.099
refrel = 0.75
ang = 0.0
nang = 7

# ----- Python -----
S1, S2, Qext, Qsca, Qback, gsca = bhmie_herbert_kaiser_july2012.bhmie(x, refrel, array([ang]))
print("bhmie, Python implementation:")
for i in range(len(S1)):
    print("{}  {}".format(S1, S2))

# ----- FORTRAN -----
S1, S2, Qext, Qsca, Qback, Gsca = bhmie_fortran.bhmie(x, refrel, nang)
print("bhmie, FORTRAN implementation:")
for i in range(2 * nang - 1):
    print("{}  {}  ".format(S1[i], S2[i]))
print(Qext)
print(Qback)
print(Gsca)

# ----- MATLAB, one file -----
print('Calling MATLAB')
S1, S2 = call_bhmie_matlab(x, refrel, nang)
print('bhmie, MATLAB single-file implementation')
print(S1, S2)

# ----- MATLAB, multi-file -----
print('Calling MATLAB')
S1, S2 = call_Mie_MATLAB(refrel, x, ang)
print('Multi-file MATLAB implementation')
print(S1, S2)
