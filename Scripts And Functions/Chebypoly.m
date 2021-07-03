function psi = Chebypoly(variable,n,a,b)
%% Inputs:
% variable - wanted coordinate
% n - number of Chebyshev polynomials to use
% a - lower region bound
% b - upper region bound
%% Example:
% syms x real, n=4; a=3; b=13; 
% psy=Chebypoly(x,n,a,b).';
% x1=linspace(a,b,100)'; % define 100 points from a to b
% fpsy = matlabFunction(psy+eps*x); % make a matlab function
% plot(x1, fpsy(x1))
% shg

x=sym(variable);
psi=chebyshevT(0:n-1,(2*(x-a)-(b-a))/(b-a));
psi=psi';
end
