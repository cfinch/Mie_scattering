import os
from numpy import *
import h5py
from bhmie_fortran import bhmie_fortran
from Mie_MATLAB.Mie_MATLAB_interface import call_Mie_MATLAB
from bhmie_matlab.bhmie_matlab_interface import call_bhmie_matlab
import bhmie_herbert_kaiser_july2012

def prettyprint(S1, S2, Qext, Qsca, Qback, gsca, nang):
    """Display data in a layout similar Wiscombe's format"""
    print("Angle  Cosine              S1                        S2             Intensity         Polzn")
    for i in range(nang):
        print("{:.2f}  {: .4f}  {: .5e} {: .5e}  {: .5e} {: .5e}     TBD             TBD".format(angles[i],
            cos(angles[i]), S1[i].real, S1[i].imag, S2[i].real, S2[i].imag))
        print("                ({: .5f})   ({: .5f})    ({: .5f})   ({: .5f})".format(
            S1[i].real / TestS1[case, i].real, S1[i].imag / TestS1[case, i].imag,
            S2[i].real / TestS2[case, i].real, S2[i].imag / TestS2[case, i].imag))

# Setup
fileName = "Wiscombe/Wiscombe_data.h5"

# ----- Read Wiscombe's validation data -----
inFile = h5py.File(fileName, 'r')

# S1 and S2, complex
TestS1 = empty(inFile['S1'].shape, dtype=complex128)
TestS2 = empty(inFile['S2'].shape, dtype=complex128)
for i in range(inFile['S1'].shape[0]):
    for j in range(inFile['S1'].shape[1]):
        TestS1[i,j] = inFile['S1'][i,j][0] + 1j * inFile['S1'][i,j][1]
        TestS2[i,j] = inFile['S2'][i,j][0] + 1j * inFile['S2'][i,j][1]

# refractive index and forward/backward scattering, complex
TestCR = array( [complex128(num[0] + 1j * num[1]) for num in inFile['TestCR']], dtype=complex128)
TestSF = array( [complex128(num[0] + 1j * num[1]) for num in inFile['TestSF']], dtype=complex128)
TestSB = array( [complex128(num[0] + 1j * num[1]) for num in inFile['TestSB']], dtype=complex128)

# Real arrays
TestXX = inFile['TestXX'][:]    # size parameter
TestQE = inFile['TestQE'][:]    # QEXT 
TestQS = inFile['TestQS'][:]    # QSCA
TestGQ = inFile['TestGQ'][:]    # GQSC

inFile.close()


angles = arange(0.0, pi + 0.01, pi/6)
nang = len(angles)

# Loop over Wiscombe's test cases
for case in range(4, len(TestCR)):
    refrel = TestCR[case]
    x = TestXX[case]
    
    print("---------------------------------------------------------------" + os.linesep)
    print("Case {}".format(case + 1))   # match Wiscombe's case numbers
    print("Refractive index = {:.4f}".format(refrel))
    print("Size parameter = {:.4f}".format(x) + os.linesep)

    # ----- Python -----
    S1, S2, Qext, Qsca, Qback, gsca = bhmie_herbert_kaiser_july2012.bhmie(x, refrel, angles)
    print("bhmie, Python implementation:")
    prettyprint(S1, S2, Qext, Qsca, Qback, gsca, len(angles))

    # ----- FORTRAN -----
    S1, S2, Qext, Qsca, Qback, gsca = bhmie_fortran.bhmie(x, refrel, 4)
    print("bhmie, FORTRAN implementation:")
    prettyprint(S1, S2, Qext, Qsca, Qback, gsca, nang)

    # ----- MATLAB, one file -----
    print('Calling MATLAB single-file implementation')
    for a in range(len(angles)):
        S1[a], S2[a] = call_bhmie_matlab(x, refrel, array(angles[a]))
    prettyprint(S1, S2, 0.0, 0.0, 0.0, 0.0, nang)
    
    # ----- MATLAB, multi-file -----
    print('Calling MATLAB multi-file')
    for a in range(len(angles)):
        S1[a], S2[a] = call_Mie_MATLAB(x, refrel, angles[a])
    prettyprint(S1, S2, 0.0, 0.0, 0.0, 0.0, nang)

