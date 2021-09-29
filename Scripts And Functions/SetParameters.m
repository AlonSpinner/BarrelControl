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
npsi=6; %number of base functions
nkill4kalman=2; %number of modes to kill for kalman (zerofy pertubations)
ts_controller=1/1000; %[s]

minModalDamp=2*zeta*(30*2*pi)/1e13;

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
[Gama,Lambda,PHI]=CalcualteMatrices(psi,dpsi,ddpsi,L,E,I,A,rho,zeta,npsi,C,minModalDamp);

g=9.81;
Fg=double(int(-g*rho*A*psi,xi,0,L));

F2Q_xi0=double(subs(psi,xi,0));
M2Q_xi0=double(subs(dpsi,xi,0));
F2Q_xiL=double(subs(psi,xi,L));

psi_0=double(subs(psi,0));
dpsi_0=double(subs(dpsi,0));
psi_Ldiv8=double(subs(psi,L/8));
dpsi_Ldiv8=double(subs(dpsi,L/8));
psi_L=double(subs(psi,L));
dpsi_L=double(subs(dpsi,L));
%% Calculating Kalman Matrices
%first mode is free mode - wr=0
%x = [eta(2:end), etadot(2:end)]'
%xdot= [etadot(2:end), etaddot(2:end)]'

%Continous system:
%measurement y_hat - theta_dot(xi=0,t) where theta=w', w is displacement
E=1*E; %lie for kalman filter! :)
A=1*A; %lie for kalman filter! :)
[kGama,kLambda,kPHI]=CalcualteMatrices(psi,dpsi,ddpsi,L,E,I,A,rho,zeta,npsi,C,minModalDamp);

nk=npsi-nc;

kalman_A=[zeros(nk), eye(nk);
        -kLambda, -kGama];
kalman_B=[zeros(nk);eye(nk)];
kalman_C=[psi_0'*C*kPHI*[eye(nk),zeros(nk)];
    dpsi_Ldiv8'*C*kPHI*[eye(nk),zeros(nk)]];
kalman_D=zeros(size(kalman_C,1),nk);
kalman_sys=ss(kalman_A,kalman_B,kalman_C,kalman_D);

%minimal realization and make it discrete
d_kalman_sys=c2d(kalman_sys,ts_controller);

kmr_A=d_kalman_sys.A;
kmr_B=d_kalman_sys.B;
kmr_C=d_kalman_sys.C;
kmr_D=d_kalman_sys.D;

Nw=size(kmr_A,2); %amount of values in state
Ny=size(kmr_C,1); %measurements size
kmr_Q = eye(Nw,Nw); %process covaraince matrix
kmr_R = eye(Ny); %measurement covariance matrix
kmr_N = zeros(Nw,Ny);%noise cross covaraince matrix, between process and measurement noise

eta0=zeros(Nw,1);
%% Assign Paramters to the Model
mdlWks = get_param('BarrelControl','ModelWorkspace');
clear(mdlWks)

%'Raw' parameters
assignin(mdlWks,'L',L);

%Dynamics
assignin(mdlWks,'ts_controller',ts_controller); 
assignin(mdlWks,'C',C);
assignin(mdlWks,'PHI',PHI);
assignin(mdlWks,'Lambda',Lambda);
assignin(mdlWks,'Gama',Gama);

assignin(mdlWks,'psi_0',psi_0);
assignin(mdlWks,'dpsi_0',dpsi_0);
assignin(mdlWks,'dpsi_Ldiv8',dpsi_Ldiv8);
assignin(mdlWks,'psi_L',psi_L);
assignin(mdlWks,'dpsi_L',dpsi_L);

assignin(mdlWks,'Fg',Fg);
assignin(mdlWks,'F2Q_xi0',F2Q_xi0);
assignin(mdlWks,'F2Q_xiL',F2Q_xiL);
assignin(mdlWks,'M2Q_xi0',M2Q_xi0);

assignin(mdlWks,'kmr_A',kmr_A);
assignin(mdlWks,'kmr_B',kmr_B);
assignin(mdlWks,'kmr_C',kmr_C);
assignin(mdlWks,'kmr_D',kmr_D);
assignin(mdlWks,'kmr_Q',kmr_Q);
assignin(mdlWks,'kmr_R',kmr_R);
assignin(mdlWks,'kmr_N',kmr_N);
assignin(mdlWks,'eta0',eta0);

assignin(mdlWks,'npsi',npsi);
assignin(mdlWks,'nc',nc);
assignin(mdlWks,'nk',nk);
assignin(mdlWks,'nkill4kalman',nkill4kalman);
%% need to write down manualy:
proj = matlab.project.rootProject;
rootFolder=proj.RootFolder;
f_psi=matlabFunction(psi,'File',fullfile(rootFolder,'Scripts And Functions','computePsi')); %to drawing Sfcn

function [Gama,Lambda,s_PHI]=CalcualteMatrices(psi,dpsi,ddpsi,L,E,I,A,rho,zeta,npsi,C,minModalDamp)
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

%minModalDamp is important as without it we have a pure generalized forces integration mode.
%which matlab can't handle in the kalman filter.
%We could either do this or change the shape function prehaps? this seems
%better
damping=max(2*zeta*s_wr,minModalDamp); 
Gama=diag(damping);
Lambda=diag(s_wr.^2);
end