%In order to calculate the Magnetic Dipole we need to first calculate the
%demagnetization factor Nd


function [N_d] = Demag_factor(L, r)
% This function calculates the value of N_d according to the formula:
% N_d = 4[ln(L/r) - 1] / ((L/r)^2 - 4*ln(L/r))

% Input:
% L: The value of L.
% r: The value of r.

% Output:
% N_d: The calculated value of N_d.

N_d = 4*(log(L/r) - 1) / ((L/r)^2 - 4*log(L/r));
end


