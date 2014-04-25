% GNR Model
%close all
clear all

Temp_K = 300; % (K) Temperature

num_layers = 8; % number of layers of graphene in GNR

gnr_length = 20e-6; % (m) length of GNR XC
gnr_width = 7.8e-9; % (m) width of GNR XC
h_vert = 1.99e-6;

rho_interlayer = 3e-1; %(Ohm*m) c-axis resistivity of GNR - Between 3e-3 and 3e-1 Ohm*m
Ef = 0.2; % (eV) Fermi energy of the bottom GNR layer
d_m = 0.35e-9; % (m) interlayer separation in GNR

Rq = 12.9e3; % (Ohm) Quantum resistance of a single-layer GNR
Rc = 0; % (Ohm) Contact resustance between contact pads and GNR XC
vf = 8e5; % (m/s) Fermi velocity in graphene
mfp_defect = 1000e-9; % (m) MFP from defect scattering in graphene on SiO2 (100-300nm realistic, 1um ideal)
prob_backscatter = 0.0; % (-) Backscattering probability ( 0 for ballistic, 0.2 realistic)
epsrd = 16.5; % ILD relative permittivity

kb = 1.381e-23; % Boltzmann constant (J/K)
h0 = 6.626e-34; % (J*s) Planck constant
q = 1.602e-19; % (C) electron charge 
h = h0/q; % (eV*s) Planck constant
eps0 = 8.854e-12; % (F/m) Vacuum permittivity

kbT = kb/q*Temp_K; % (eV) combined Boltzmann temp factor
epsd = epsrd*eps0;


%% Calculate number of conducting channels
% beta = 0; % 0 for metallic, 1/3 for semiconducting
% 
% Esubm = @(m) h*vf/2/W * abs(m+beta);
% Nchc = @(m) (1 + exp((Esubm(m) - Ef)/kbT) )^-1;
% Nchv = @(m) (1 + exp((Esubm(m) + Ef)/kbT) )^-1;
% Nchm = @(m) Nchc(m) + Nchv(m);
% 
% keep_going = 1;
% Nch = 0;
% m = 1;
% 
% while(keep_going)
%     Nch = Nch + Nchm(m);
%     keep_going = (Nchm(m) > 0.01);
%     m = m+1;
% end
% m_max = m-1;

gnr_width_nm = gnr_width*1e9;
mfp_defect_nm = mfp_defect * 1e9;
[ Nch,Rsheet,leff] = sheetres_single_mGNR_ld1u_mod( Ef,gnr_width_nm,prob_backscatter,Temp_K,mfp_defect_nm);   

%% Calculate effective mean free path
% 
% lambda_edge_m = @(E,m) W/P*sqrt( (E/Esubm(m))^2-1);
% f = @(E) (1 + exp((E-Ef)/kbT)).^-1;
% 
% lambda_eff_vec = zeros(1,m_max);
% lambda_edge_vec = zeros(1,m_max);
% for m = 1:m_max
%     lambda_eff_vec(m) = (1/lambda_d + 1/lambda_edge_m(Ef,m))^-1;
%     lambda_edge_vec(m) = lambda_edge_m(Ef,m);
% end
% 
% num = 0;
% den = 0;
% for m = 1:m_max
%     num = num + lambda_eff_vec(m)*f(Esubm(m));
%     den = den + f(Esubm(m));
% end
% 
% lambda_eff_calc = num/den;

%% Nevermind, just use an effective value

% lambda_eff = lambda_d;
lambda_eff = leff*1e-9;


%% Calculate component resistances

ry = rho_interlayer*d_m/gnr_width/gnr_length;
rx = Rq/Nch/lambda_eff*gnr_length;

%% Top contacted MLGNR resistance

rn = zeros(1,num_layers);
rn(1) = rx;
for n = 2:num_layers
    ra = rn(n-1) + ry;
    rb = rx;
    rn(n) = (ra^-1 + rb^-1)^-1;
end

r_tc = rn(end);

%% Side-contacted MLGNR resistance

r_sc = rx/num_layers;

%% Inductance

Lk = h0/4/q^2/vf/num_layers/Nch; %(H/m) Kinetic inductance
L = Lk;

%% Capacitance

M1 = @(a) 2*pi./( log( 2*(1+(1-a.^2).^(1/4))./(1-(1-a.^2).^(1/4)) ));
M2 = @(a) 2/pi*log( (2*(1+sqrt(a)))./(1-sqrt(a)) );

M = @(a) ((0 <= a) & (a<= 1/sqrt(2))) .* M1(a) + ...
         ((1/sqrt(2) < a) & (a<= 1)) .* M2(a);

Cq = 4*q^2/h0/vf * num_layers * Nch; % (F/m) Quantum capacitance;
Ce = epsd*M(tanh(pi*gnr_width/4/h_vert));

c = (Cq^-1 + Ce^-1)^-1;

%% Delay

Rtc = r_tc ;
Rsc = r_sc ;

rtc_pul = Rtc/gnr_length;
rsc_pul = Rsc/gnr_length;

C = c * gnr_length;
L = Lk * gnr_length;

wn = 1/gnr_length*sqrt(2/Lk/c);
alpha_tc = rtc_pul/2/wn/Lk;
alpha_sc = rsc_pul/2/wn/Lk;

Rtc
Rsc

delay = @(wn,gamma) (1.047*exp(gamma/0.85) + 1.39*gamma)/wn;
delay_sc = delay(wn,alpha_sc)
delay_tc = delay(wn,alpha_tc)


