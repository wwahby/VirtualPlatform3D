% GNR Model
close all
clear all

T = 300; % (K) Temperature

M = 100; % number of segments to break XC into
N = 8; % number of layers of graphene in GNR

length = 20e-6; % (m) length of GNR XC
W = 10e-9; % (m) width of GNR XC
h_vert = 1.99e-6;

rho_c = 3e-3; %(Ohm*m) c-axis resistivity of GNR - Between 3e-3 and 3e-1 Ohm*m
Ef = 0.2; % (eV) Fermi energy of the bottom GNR layer
d_m = 0.35e-9; % (m) interlayer separation in GNR

Rq = 12.9e3; % (Ohm) Quantum resistance of a single-layer GNR
Rc = 0; % (Ohm) Contact resustance between contact pads and GNR XC
vf = 8e5; % (m/s) Fermi velocity in graphene
lambda_d = 1000e-9; % (m) MFP from defect scattering in graphene on SiO2 (100-300nm realistic, 1um ideal)
P = 0.2; % (-) Backscattering probability ( 0 for ballistic, 0.2 realistic)
epsrd = 16.5; % ILD relative permittivity

kb = 1.381e-23; % Boltzmann constant (J/K)
h0 = 6.626e-34; % (J*s) Planck constant
q = 1.602e-19; % (C) electron charge 
h = h0/q; % (eV*s) Planck constant
eps0 = 8.854e-12; % (F/m) Vacuum permittivity

kbT = kb/q*T; % (eV) combined Boltzmann temp factor
dx = length/M;
epsd = epsrd*eps0;


%% Calculate number of conducting channels
beta = 0; % 0 for metallic, 1/3 for semiconducting

Esubm = @(m) h*vf/2/W * abs(m+beta);
Nchc = @(m) (1 + exp((Esubm(m) - Ef)/kbT) )^-1;
Nchv = @(m) (1 + exp((Esubm(m) + Ef)/kbT) )^-1;
Nchm = @(m) Nchc(m) + Nchv(m);

keep_going = 1;
Nch = 0;
m = 1;

while(keep_going)
    Nch = Nch + Nchm(m);
    keep_going = (Nchm(m) > 0.01);
    m = m+1;
end
m_max = m-1;

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

lambda_eff_calc = num/den;

%% Nevermind, just use an effective value

lambda_eff = lambda_d;


%% Calculate component resistances

ry = rho_c*d_m/W/length;
rx = Rq/Nch/lambda_eff*length;

%% Top contacted MLGNR resistance

rn = zeros(1,N);
rn(1) = rx;
for n = 2:N
    ra = rn(n-1) + ry;
    rb = rx;
    rn(n) = (ra^-1 + rb^-1)^-1;
end

r_tc = rn(end);

%% Side-contacted MLGNR resistance

r_sc = rx/N;

%% Inductance

Lk = h0/4/q^2/vf/N/Nch; %(H/m) Kinetic inductance
L = Lk;

%% Capacitance

M1 = @(a) 2*pi*( log( 2*(1+(1-a^2)^(1/4))/(1-(1-a^2)^(1/4))))^-1;
M2 = @(a) 2/pi*log( (2*(1+sqrt(a)))/(1-sqrt(a)) );

M = @(a) ((0 <= a) && (a<= 1/sqrt(2))) * M1(a) + ...
         ((1/sqrt(2) < a) && (a<= 1)) * M2(a);

Cq = 4*q^2/h0/vf * N * Nch; % (F/m) Quantum capacitance;
Ce = epsd*M(tanh(pi*W/4/h_vert));

c = (Cq^-1 + Ce^-1)^-1;

%% Delay

Rtc = r_tc ;
Rsc = r_sc ;

C = c * length;
L = Lk * length;

wn = 1/length*sqrt(2/L/C);
alpha_tc = Rtc/2/wn/L;
alpha_sc = Rsc/2/wn/L;

Rtc
Rsc

tau_sc = (1.047*exp(alpha_sc/0.85) + 1.39*alpha_sc)/wn
tau_tc = (1.047*exp(alpha_tc/0.85) + 1.39*alpha_tc)/wn


