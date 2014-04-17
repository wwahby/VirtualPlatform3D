% GNR Model - Kumar 2012
close all
clear all

T = 300; % (K) Temperature

M = 1e1; % number of segments to break XC into
N = 6; % number of layers of graphene in GNR

L = 20e-6; % (m) length of GNR XC
W = 10e-9; % (m) width of GNR XC

rho_c = 3e-2; %(Ohm*m) c-axis resistivity of GNR - Between 3e-3 and 3e-1 Ohm*m
Ef = 0.2; % (eV) Fermi energy of the bottom GNR layer
d_m = 0.35e-9; % (m) interlayer separation in GNR

Rq = 12.9e3; % (Ohm) Quantum resistance of a single-layer GNR
Rc = 0; % (Ohm) Contact resustance between contact pads and GNR XC
vf = 8e5; % (m/s) Fermi velocity in graphene
lambda_d = 1000e-9; % (m) MFP from defect scattering in graphene on SiO2
P = 0.0; % (-) Backscattering probability

kb = 1.381e-23; % Boltzmann constant (J/K)
h0 = 6.626e-34; % (J*s) Planck constant
q = 1.602e-19; % (C) electron charge 
h = h0/q; % (eV*s) Planck constant

kbT = kb/q*T; % (eV) combined Boltzmann temp factor
dx = L/M;


%% Calculate number of conducting channels
beta = 0;

Esubm = @(m) h*vf/2/W * abs(m+beta);
Nchm = @(m) (1 + exp((Esubm(m) - Ef)/kbT) )^-1;

keep_going = 1;
Nch = 0;
m = 1;
% while(keep_going)
%     if (Esubm(m) - Ef < 0)
%         Nch = Nch + Nchm(m);
%     
%         keep_going = (Nchm(m) > 0.05);
%         m = m+1;
%     else
%         keep_going = 0;
%     end
% end

while(keep_going)
    Nch = Nch + Nchm(m);
    keep_going = (Nchm(m) > 0.01);
    m = m+1;
end
m_max = m-1;

W_nm = W*1e9;
[ Nch,Rsheet,leff] = sheetres_single_mGNR_ld1u( Ef,W_nm,P,T);
%% Calculate effective mean free path

lambda_edge_m = @(E,m) W/P*sqrt( (E/Esubm(m))^2-1);
f = @(E) (1 + exp((E-Ef)/kbT)).^-1;

lambda_eff_vec = zeros(1,m_max);
lambda_edge_vec = zeros(1,m_max);
for m = 1:m_max
    lambda_eff_vec(m) = (1/lambda_d + 1/lambda_edge_m(Ef,m))^-1;
    lambda_edge_vec(m) = lambda_edge_m(Ef,m);
end

num = 0;
den = 0;
for m = 1:m_max
    num = num + lambda_eff_vec(m)*f(Esubm(m));
    den = den + f(Esubm(m));
end

lambda_eff = num/den;
lambda_eff = lambda_d;


%% Calculate component resistances
Ra = rho_c*d_m/W/dx;
Rb = Rq*dx/Nch/lambda_eff;

%% Construct H matrices
h11d0 = 1+2*Ra/Rb*ones(1,M);
h11dp1 = -Ra/Rb*ones(1,M);
h11dm1 = -Ra/Rb*ones(1,M);

h12d0 = 2*Ra*ones(1,M);
h12dp1 = -Ra*ones(1,M);
h12dm1 = -Ra*ones(1,M);

h21d0 = 1/Rb*ones(1,M);

h22d0 = 1*ones(1,M);

H11 = spdiags( [h11dm1' h11d0' h11dp1'], [-1 0 1], M, M );
H12 = spdiags( [h12dm1' h12d0' h12dp1'], [-1 0 1], M, M );
H21 = spdiags( h21d0' , 0, M, M );
H22 = spdiags( h22d0' , 0, M, M );


H = [H11 H12 ; H21 H22];
A = H^(N-1);

%% Pull the sub-matrices out
A11 = A(1:M, 1:M);
A12 = A(1:M, M+1:2*M);
A21 = A(M+1:2*M, 1:M);
A22 = A(M+1:2*M, M+1:2*M);

%% calculate effective resistance
%Reff_mat =  A11 * (1/Rb*A11 - A21)^-1;
Reff_mat =  A11 * (1/Rb*A11 + A21)^-1;

Reff = full(sum(sum(Reff_mat)))


%Rnet = Reff + (Rq + 2*Rc)/Nch


