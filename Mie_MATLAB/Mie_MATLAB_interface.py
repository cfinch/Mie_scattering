import subprocess
import os
import ast
import math

def call_Mie_MATLAB(size_parameter, refractive_index, angle):
    """Usage: S1, S2 = call_Mie_MATLAB(0.75, 0.101, 3.14159/6)"""

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
    matlab_script += "S12 = Mie_S12({}, {}, {}); ".format(refractive_index, 
        size_parameter, math.cos(angle))
    matlab_script += "disp(['S1 = ', num2str(real(S12(1))), ' ', num2str(imag(S12(1)))]); "
    matlab_script += "disp(['S2 = ', num2str(real(S12(2))), ' ', num2str(imag(S12(2)))]); "
    matlab_script += "exit"

    process = subprocess.Popen(["/opt/MATLAB/R2011a/bin/matlab", "-nojvm",
        "-nodisplay", "-r", matlab_script],stdout=subprocess.PIPE,
        cwd="/home/cfinch/Projects/Mie_scattering/Mie_MATLAB")

    # Get and parse text output from MATLAB
    output = process.stdout.readlines()
    for line in output:
        if line.find('S1') > -1:
            S1 = parse_output(line)
        if line.find('S2') > -1:
            S2 = parse_output(line)

    return S1, S2
