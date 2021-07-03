%% Initalize 
syms xi real
L=5; %m
rho=8e3; %kg/m^3
E=1*210e9; %210Gpa
rm_in=15e-3;
rm_out=25e-3- 5e-3*xi/L; %m, symbolic

%Calculating area and moment of inertia
A=pi*(rm_out^2 -rm_in^2); %m^2
I=pi/4*(rm_out^4-rm_in^4); %m^2

Mass=double(int(A*rho,0,L));
%% Calculating Matrices
%Number of base functions for rayligh ritz
n=8;

psi=Chebypoly(xi,n,0,L);
dpsi=diff(psi,xi);
ddpsi=diff(dpsi,xi);

K=double(E*int(I*(ddpsi*ddpsi'),0,L));
M=double(rho*int(A*(psi*psi'),0,L));

c = double([
%     [subs(psi',xi,0)]; %w(0)=0
    [subs(psi',xi,L/8)]; %w(L/8)=0
%     [subs(dpsi',xi,0)]; %w'(0)=0
            ]);
C=null(c,'r');

%Apply constraints on matrices
Kc=C'*K*C;
Mc=C'*M*C;

%ensure symmatrey
Kc=0.5*(Kc+Kc');
Mc=0.5*(Mc+Mc');

[PHI,wr2]=eig(Kc,Mc); %solves Av=lambda*Bv ~ K-(wn^2)*M
wr=sqrt(diag(wr2));
wrMat(1:length(wr),n)=wr;

[s_wr,Indx]=sort(wr); %sort
s_wr=abs(s_wr'); %abs or real?
s_PHI=PHI(:,Indx); %PHI to work with

zeta=0.03;
Gama=diag(2*zeta*s_wr);
Lambda=diag(s_wr.^2);

g=9.81;
Fg=double(int(-g*psi,xi,0,L));

F2Q_xi0=double(subs(psi,xi,0));
M2Q_xi0=double(subs(dpsi,xi,0));
%% Assign Paramters to the Model
mdlWks = get_param('BarrelControl','ModelWorkspace');
clear(mdlWks)

%'Raw' parameters
assignin(mdlWks,'L',L);

%Dynamics
assignin(mdlWks,'C',C);
assignin(mdlWks,'PHI',s_PHI);
assignin(mdlWks,'Lambda',Lambda);
assignin(mdlWks,'Gama',Gama);
assignin(mdlWks,'psi_0',double(subs(psi,0)));
assignin(mdlWks,'dpsi_0',double(subs(dpsi,0)));
assignin(mdlWks,'psi_L',double(subs(psi,L)));
assignin(mdlWks,'dpsi_L',double(subs(dpsi,L)));
assignin(mdlWks,'Fg',Fg);
assignin(mdlWks,'F2Q_xi0',F2Q_xi0);
assignin(mdlWks,'M2Q_xi0',M2Q_xi0);
%% need to write down manualy:
f_psi=matlabFunction(psi); %to drawing Sfcn