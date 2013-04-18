% MIEV0 Test Case  6 from NCAR Technical Note
% Mie Scattering Calculations: Advances in Technique and Fast, Vector-Speed Computer Codes
% Warren J. Wiscombe, 1996
% For unknown reasons, Wiscombe's code returns the complex conjugate of the
% S1 and S2 values reported by codes based on Bohren and Huffman

angles = [0.0:pi/6:pi];
x = 0.1010;
m = 0.75;

disp(' Angle                 S1                             S2');
for i=1:length(angles)
    S12 = Mie_S12(m, x, cos(angles(i)));
    str = sprintf('%6.2f    % e % e    % e % e', angles(i)*180/pi, real(S12(1)), imag(S12(1)), real(S12(2)), imag(S12(2)));
    disp(str);
end
