function psi = computePsi(xi)
%COMPUTEPSI
%    PSI = COMPUTEPSI(XI)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    24-Sep-2021 09:17:00

t2 = xi.^2;
psi = [1.0;xi.*(2.0./3.0)-1.0;t2.*(8.0./9.0)-xi.*(8.0./3.0)+1.0;t2.*(-1.6e+1./3.0)+xi.*6.0+xi.^3.*(3.2e+1./2.7e+1)-1.0];
