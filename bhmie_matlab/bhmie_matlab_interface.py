import subprocess
import os
import ast

def call_bhmie_matlab(size_parameter, refractive_index, angle):
    """Usage: S1, S2 = call_bhmie_matlab(0.75, 0.101, pi/2)"""

    def parse_output(string):
        """Parse text output from MATLAB with the following content:
        S2 = -1.8176e-08 -0.00016481
        Returns a Python complex number
        """
        words = line.split()
        real_string = words[-2]
        imag_string= words[-1]
        real = ast.literal_eval(real_string)
        imag = ast.literal_eval(imag_string)

        number = complex(real, imag)

        return number

    # Construct MATLAB script
    matlab_script = "disp(['running MATLAB']); "
    matlab_script += "[s1, s2, qext, qsca, qback, qgsca] = mie({}, {}, {}); ".format(
        size_parameter, refractive_index, angle)
    matlab_script += "disp(['S1 = ', num2str(real(s1)), ' ', num2str(imag(s1))]); "
    matlab_script += "disp(['S2 = ', num2str(real(s2)), ' ', num2str(imag(s2))]); "
    matlab_script += "exit"

    process = subprocess.Popen(["/opt/MATLAB/R2011a/bin/matlab", "-nojvm",
        "-nodisplay", "-r", matlab_script],stdout=subprocess.PIPE,
        cwd="/home/cfinch/Projects/Mie_scattering/bhmie_matlab")

    # Get and parse text output from MATLAB
    output = process.stdout.readlines()
    for line in output:
        if line.find('S1') > -1:
            S1 = parse_output(line)
        if line.find('S2') > -1:
            S2 = parse_output(line)

    return S1, S2
