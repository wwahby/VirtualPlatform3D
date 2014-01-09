% Power and signal codesign

%% Chip descriptors
% Stack parameters
Ng = 1.0e9;
S = 4;
Ach_mm2 = 100;
Ach_m2 = Ach_mm2*1e-6;

% gate parameters
eps_ox = 25; % HfO2
tox = 1e-9;
w_trans = 32e-9;
N_trans_per_gate = 4;
Ioff = 10e-9; %(A/um)

% Tsv parameters
Atf_max = 0.01; % maximum allowable TSV area, as a fraction of total chip area
gate_pitch = 300e-9;
AR_tsv = 20;

% Rent parameters
p = 0.6; % rent exponent
fo = 4; % avg fanout
alpha = fo/(fo+1); % input terminal fraction
k = 3/alpha; %rent constant

% Wiring parameters
fmax = 1e9;
chi = 2/3;
rho_m = 17.2e-9; % Cu
epsr_d = 3.0; % Low-k dielectric
Tclk = 1/fmax; % (s)
alpha_t = 1.1*6.2;

% Repeater parameters
Ro = 1e3; % (Ohm) Gate output resistance

% Power parameters
a = 0.1; % logic activity factor
Vdd = 1.0; % (V)

% Model flags
use_joyner = 0; % 0=Use corrected distribution, 1=use original joyner distribution
redo_wiring = 0; % Redo wire layer assignment after repeater insertion? 0=no, 1=yes

%% Power supply noise model parameters
psn_target = 0.075;   % (mV) Acceptable power supply noise
mu_m = 1.257e-6;      %copper permeability

decap = 0.1; % (Ratio) - Fraction of chip area dedicated to decoupling capacitors
Npad1d = 32;
Npad = Npad1d^2;         %total number of power or ground pads
Ngrid = 21*21;        %grid fineness
padsize = 1;          % Pad size is in terms of segment number
Tseg = 1e-6;          % Thickness of grid segment
Wseg = 2e-6;          % width of grid segment

T = 25; % (deg C) Temperature

% Package parameters
RPKG = 0.006; % (Ohm)
LPKG = 0.5e-9; % (H)

tic % begin timing
%% TSV number determination
% Inputs:
%   Rent parameters
%   Number of logic gates
%   Number of layers
disp(' ')
disp('Determining system parameters...')
[nt_max nt_tot nt_to nt_through Tacmat] = estimate_tsvs_required(Ng,S,k,p,alpha);

%% TSV Sizing
% Inputs:
%   Area available
%   Max area for TSVs
%   TSV aspect ratio
Ach = Ng/S;
tsv_area_ratio = Atf_max;
[w_tsv_gp h_tsv_gp] = size_tsvs(Ach, tsv_area_ratio, nt_max, AR_tsv );
%h_tsv_gp = round(h_tsv_gp);
h_tsv_m = h_tsv_gp*gate_pitch;
w_tsv_m = w_tsv_gp*gate_pitch;

%% System determination
% Run WLD + WLA + RI to get power estimate
[ iidf l Ln pn pn_orig Cxc Ltot Cn Pdyn Plk Pw Prep Ng_act N_tsvs T_tsvs Atf_act ] = gen_design(Ng,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);
Ptot = Pdyn + Plk + Pw + Prep;
Pdens = Ptot/Ach_m2;
pitch_tsv = T_tsvs*gate_pitch*100; % [FIX] Power TSVs aren't going to be on the same pitch as signal TSVs

%% Power noise estimation
% Inputs:
%   Total power (from previous step)
%   TSV geometry from TSV Sizing step
disp('Evaluating power supply noise...')

psn_iterations = 1;
psn_max = calc_psn(Ach_m2,Npad,rho_m,mu_m,S,Vdd,Pdens,decap,h_tsv_m,w_tsv_m,pitch_tsv,RPKG,LPKG,Ngrid,padsize,Tseg,Wseg,T);
dispstr = sprintf('\tpsn_runs: %d\tNpads: %d\tpsn_target: %d\tpsn_max: %d',psn_iterations,Npad, psn_target, psn_max);
disp(dispstr)

% These three while loops help us hone in on the actual psn target
% This is a pretty ghetto way of doing it, but it's quick and gets us a lot
% closer than just doing one loop (and was quicker to code than a binary
% search or what have you)
%   The idea is to just blindly increase the number of power/gnd pads by some
% multiplicative factor until we satisfy the power noise constraint
%   Then we decrease the power/gnd pad allocation by a (smaller) factor
%  until we no longer satisfy the constraint
%   And finally we increase the pad allocation again (by an even smaller
%   factor) until we satisfy the constraint again
% This method is definitely slower than something more intelligent (and probably doesn't get us as close to to the target as possible), but
% since the PSN calculation runs reasonably fast and since we can get close
% enough to the PSN target (and since the power/gnd tsv demands are usually
% much lower than the signal TSV requirements), we can live with these
% inefficiencies for now

%Overshoot
while ((psn_max > psn_target) && (psn_iterations < 20))
    Npad1d = round(1.3*Npad1d);
    Npad = Npad1d^2;
    psn_max = calc_psn(Ach_m2,Npad,rho_m,mu_m,S,Vdd,Pdens,decap,h_tsv_m,w_tsv_m,pitch_tsv,RPKG,LPKG,Ngrid,padsize,Tseg,Wseg,T);
    psn_iterations = psn_iterations+1;
    
    dispstr = sprintf('\tpsn_runs: %d\tNpads: %d\tpsn_target: %d\tpsn_max: %d',psn_iterations,Npad, psn_target, psn_max);
    disp(dispstr)
end

%Undershoot
while ((psn_max < psn_target) && (psn_iterations < 20))
    Npad1d = round(0.9*Npad1d);
    Npad = Npad1d^2;
    psn_max = calc_psn(Ach_m2,Npad,rho_m,mu_m,S,Vdd,Pdens,decap,h_tsv_m,w_tsv_m,pitch_tsv,RPKG,LPKG,Ngrid,padsize,Tseg,Wseg,T);
    psn_iterations = psn_iterations+1;
    
    dispstr = sprintf('\tpsn_runs: %d\tNpads: %d\tpsn_target: %d\tpsn_max: %d',psn_iterations,Npad, psn_target, psn_max);
    disp(dispstr)
end

%Overshoot
while ((psn_max > psn_target) && (psn_iterations < 20))
    Npad1d = round(1.05*Npad1d);
    Npad_tentative = Npad1d^2;
    if(Npad_tentative == Npad)
        Npad = Npad + 1;
    else
        Npad = Npad_tentative;
    end
    psn_max = calc_psn(Ach_m2,Npad,rho_m,mu_m,S,Vdd,Pdens,decap,h_tsv_m,w_tsv_m,pitch_tsv,RPKG,LPKG,Ngrid,padsize,Tseg,Wseg,T);
    psn_iterations = psn_iterations+1;
    
    dispstr = sprintf('\tpsn_runs: %d\tNpads: %d\tpsn_target: %d\tpsn_max: %d',psn_iterations,Npad, psn_target, psn_max);
    disp(dispstr)
end

% print inputs to calc_psn for debug purposes -- leave this commented out
% most of the time
% instr = sprintf('Ach_m2__%d	Npad__%d	rho_m__%d	mu_m__%d	S__%d	Vdd__%d	Pdens__%d	decap__%d	h_tsv_m__%d	w_tsv_m__%d	pitch_tsv__%d	RPKG__%d	LPKG__%d	Ngrid__%d	padsize__%d	Tseg__%d	Wseg__%d	T__%d',Ach_m2,Npad,rho_m,mu_m,S,Vdd,Pdens,decap,h_tsv_m,w_tsv_m,pitch_tsv,RPKG,LPKG,Ngrid,padsize,Tseg,Wseg,T);
% disp(instr);
%% Final report
disp(' ')
disp('Final system parameters:')
repstr = sprintf('\tNg_nom %d \t Ng_act: %d \t Atsv_nom: %.3g \t Atsv_act: %.3g \n\tN_tsvs: %d \t Npads_pow %d \t psn_nom %.4g \t psn_act %.4g', ...
                  Ng, Ng_act, Atf_max, Atf_act, N_tsvs, Npad, psn_target, psn_max);
disp(repstr)
repstr = sprintf('\th_tsv_um: %.4g \t w_tsv_um: %.4g',h_tsv_m/1e-6,w_tsv_m/1e-6);
disp(repstr);

toc % finish timing

