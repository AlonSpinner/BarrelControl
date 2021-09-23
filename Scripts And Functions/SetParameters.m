%% Initalize 
syms xi real
L=3; %m
rho=8e3; %kg/m^3
E=210e9; %210Gpa
rm_in=30e-3;
rm_out=45e-3 - 5e-3*xi/L; %m, symbolic

%Calculating area and moment of inertia
A=pi*(rm_out^2 -rm_in^2); %m^2
I=pi/4*(rm_out^4-rm_in^4); %m^2
Mass=double(int(A*rho,0,L));

%additional parameters
zeta=0.05;
npsi=8; %number of base functions

%shape function and derivative
psi=Chebypoly(xi,npsi,0,L);
dpsi=diff(psi,xi);
ddpsi=diff(dpsi,xi);

%constraints
c = double([
%     [subs(psi',xi,0)]; %w(0)=0
    [subs(psi',xi,L/8)]; %w(L/8)=0
%     [subs(dpsi',xi,0)]; %w'(0)=0
            ]);
nc=size(c,1); %number of constraints
C=null(c,'r');
%% Calculating Model Matrices
%Number of base functions for rayligh ritz
[Gama,Lambda,s_PHI]=CalcualteMatrices(psi,dpsi,ddpsi,L,E,I,A,rho,zeta,npsi,C);

g=9.81;
Fg=double(int(-g*psi,xi,0,L));

F2Q_xi0=double(subs(psi,xi,0));
M2Q_xi0=double(subs(dpsi,xi,0));
F2Q_xiL=double(subs(psi,xi,L));

%% Calculating Kalman Matrices
nk=npsi-nc;
Kalman_A=[zeros(nk), eye(nk);
        -Lambda, -Gama];
Kalman_B=blkdiag(zeros(nk),eye(nk));
Kalman_C=blkdiag(zeros(nk),eye(nk));
Kalman_D=blkdiag(zeros(nk),zeros(nk));

Ny=1; %measurements size
Nw=2*nk; %state space size
Kalman_Q = eye(Nw,Nw); %process covaraince matrix
Kalman_R = eye(Ny); %measurement covariance matrix
Kalman_N = zeros(Nw,Nw);%noise cross covaraince matrix, between process and measurement noise

eta0=zeros(Nw,1);
%% Assign Paramters to the Model
mdlWks = get_param('BarrelControl','ModelWorkspace');
clear(mdlWks)

%'Raw' parameters
assignin(mdlWks,'L',L);

%Dynamics
assignin(mdlWks,'ts_controller',1/1000); 
assignin(mdlWks,'C',C);
assignin(mdlWks,'PHI',s_PHI);
assignin(mdlWks,'Lambda',Lambda);
assignin(mdlWks,'Gama',Gama);

assignin(mdlWks,'psi_0',double(subs(psi,0)));
assignin(mdlWks,'dpsi_0',double(subs(dpsi,0)));
assignin(mdlWks,'dpsi_Ldiv8',double(subs(dpsi,L/8)));
assignin(mdlWks,'psi_L',double(subs(psi,L)));
assignin(mdlWks,'dpsi_L',double(subs(dpsi,L)));

assignin(mdlWks,'Fg',Fg);
assignin(mdlWks,'F2Q_xi0',F2Q_xi0);
assignin(mdlWks,'F2Q_xiL',F2Q_xiL);
assignin(mdlWks,'M2Q_xi0',M2Q_xi0);

assignin(mdlWks,'Kalman_A',Kalman_A);
assignin(mdlWks,'Kalman_B',Kalman_B);
assignin(mdlWks,'Kalman_C',Kalman_C);
assignin(mdlWks,'Kalman_D',Kalman_D);
assignin(mdlWks,'Kalman_Q',Kalman_Q);
assignin(mdlWks,'Kalman_R',Kalman_R);
assignin(mdlWks,'Kalman_N',Kalman_N);
assignin(mdlWks,'eta0',eta0);
assignin(mdlWks,'nk',nk);
%% need to write down manualy:
proj = matlab.project.rootProject;
rootFolder=proj.RootFolder;
f_psi=matlabFunction(psi,'File',fullfile(rootFolder,'Scripts And Functions','computePsi')); %to drawing Sfcn

function [Gama,Lambda,s_PHI]=CalcualteMatrices(psi,dpsi,ddpsi,L,E,I,A,rho,zeta,npsi,C)
K=double(E*int(I*(ddpsi*ddpsi'),0,L));
M=double(rho*int(A*(psi*psi'),0,L));

%Apply constraints on matrices
Kc=C'*K*C;
Mc=C'*M*C;

%ensure symmatrey
Kc=0.5*(Kc+Kc');
Mc=0.5*(Mc+Mc');

[PHI,wr2]=eig(Kc,Mc); %solves Av=lambda*Bv ~ K-(wn^2)*M
wr=sqrt(diag(wr2));
wrMat(1:length(wr),npsi)=wr;

[s_wr,Indx]=sort(wr); %sort
s_wr=abs(s_wr'); %abs or real?
s_PHI=PHI(:,Indx); %PHI to work with

Gama=diag(2*zeta*s_wr);
Lambda=diag(s_wr.^2);
end