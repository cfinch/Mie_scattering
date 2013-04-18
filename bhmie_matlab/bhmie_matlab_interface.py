import subprocess
import os
import ast

def call_bhmie_matlab(size_parameter, refractive_index, nang):
    """Usage: S1, S2 = call_bhmie_matlab(0.75, 0.101, 3.14159/6)"""

    def parse_output(string):
        words = line.split()
        real_string = words[-2]
        imag_string= words[-1]
        real = ast.literal_eval(real_string)
        imag = ast.literal_eval(imag_string)

        number = complex(real, imag)

        return number

    matlab_script = "disp(['running MATLAB']); [s1, s2, qext, qsca, qback, qgsca] = mie({}, {}, {});  disp(['S1 = ', num2str(real(s1)), ' ', num2str(imag(s1))]); disp(['S2 = ', num2str(real(s2)), ' ', num2str(imag(s2))]); exit".format(
        size_parameter, refractive_index, nang)

    process = subprocess.Popen(["/opt/MATLAB/R2011a/bin/matlab", "-nojvm",
        "-nodisplay", "-r", matlab_script],stdout=subprocess.PIPE,
        cwd="/home/cfinch/Projects/Mie_scattering/bhmie_matlab")

    output = process.stdout.readlines()

    for line in output:
        if line.find('S1') > -1:
#            print line
            S1 = parse_output(line)
        if line.find('S2') > -1:
#            print line
            S2 = parse_output(line)

    return S1, S2
