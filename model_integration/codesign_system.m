% Power and signal codesign
close all
clear all

%% ==================================================
%  ================== BEGIN INPUTS ==================
%  ==================================================

%% Stack parameters
S = 1;

%% 65nm Merom, entire chip
% Ng = 291e6/4;
% Ach_mm2 = 143;
% gate_pitch = 700e-9; % average gate pitch (sqrt(A_core/Ngates))
% min_pitch = 220e-9; % actual contacted gate pitch
% fmax = 3.0e9;

%% 65nm Merom, single core, wrong area
% Ng = 10e6;
% Ach_mm2 = 16;
% gate_pitch = 630e-9; % average gate pitch (sqrt(A_core/Ngates))
% min_pitch = 220e-9; % actual contacted gate pitch
% fmax = 3.0e9;

%% 65nm Merom, single core
% Ng = 10e6;
% Ach_mm2 = 27;
% gate_pitch = 820e-9; % average gate pitch (sqrt(A_core/Ngates))
% min_pitch = 220e-9; % actual contacted gate pitch
% fmax = 3.0e9;

%% 32nm Sandy Bridge, entire chip
% Ng = 2.27e9/4;
% Ach_mm2 = 435;
% gate_pitch = 875e-9;
% min_pitch = 112.5e-9;
% fmax = 3.6e9;

%% 32nm Sandy Bridge, one core
% Ng = 107e6/4;
% Ach_mm2 = 20;
% gate_pitch = 435e-9;
% min_pitch = 112.5e-9;
% fmax = 3.6e9;

%% 22nm Ivy Bridge, entire chip
% Ng = 2.86e9/4;
% Ach_mm2 = 346.5;
% gate_pitch = 348e-9;
% min_pitch = 90e-9;
% fmax = 3.0e9;

%% 22nm Ivy Bridge, one core
% Ng = 95e6/4;
% Ach_mm2 = 11.5;
% gate_pitch = 348e-9;
% min_pitch = 90e-9;
% fmax = 3.0e9;

%% Arbitrarily huge test case
Ng = 1e9;
Ach_mm2 = 100;
gate_pitch = 100e-9;
min_pitch = 100e-9;
fmax = 3.0e9;
%% Chip descriptors

Ach_m2 = Ach_mm2*1e-6;

% gate parameters
eps_ox = 25; % HfO2
tox = 1e-9;
w_trans = 30e-9;
N_trans_per_gate = 4;
Ioff = 10e-9; %(A/um)

% Tsv parameters
Atf_max = 0.10; % maximum allowable TSV area, as a fraction of total chip area
h_tsv_m_thin = 10e-6;
h_tsv_m_thick = 300e-6;
AR_tsv = 20;

% Rent parameters
p = 0.6; % rent exponent
fo = 4; % avg fanout
alpha = fo/(fo+1); % input terminal fraction
k = 3/alpha; %rent constant

% Wiring parameters
chi = 2/3;
rho_m = 17.2e-9; % Cu
epsr_d = 3.0; % Low-k dielectric
Tclk = 1/fmax; % (s)
alpha_t = 1.1*6.2;

% Repeater parameters
Ro = 1e3; % (Ohm) Gate output resistance

% Power parameters
a = 0.1; % logic activity factor
Vdd = 1.25; % (V)

%% constants
eps0 = 8.854e-12; % (F/m) vacuum permittivity

%% Power supply noise model parameters

psn.noise_target = 0.1875;              % (mV) Acceptable power supply noise
psn.decap_area_fraction = 0.1;      % (Ratio) - Fraction of chip area dedicated to decoupling capacitors
psn.Npads_1d = 32;                  % Number of pads from one side of chip to the other
psn.Npads = psn.Npads_1d^2;         % total number of power pads (does not include ground pads)
psn.Ngrid = 21*21;                  % grid fineness
psn.pad_size = 1;                   % Pad size is in terms of segment number
psn.segment_thickness = 1e-6;       % (m) Thickness of a power grid segment
psn.segment_width = 2e-6;           % (m) width of a power grid segment
psn.package_resistance = 0.006;     % (Ohm)
psn.package_inductance = 0.5e-9;    % (H)


psn.target = 0.1;   % (V) Acceptable power supply noise
mu_m = 1.257e-6;      %copper permeability

% decap = 0.1; % (Ratio) - Fraction of chip area dedicated to decoupling capacitors
% Npad1d = 32;
% Npad = Npad1d^2;         %total number of power or ground pads
% Ngrid = 21*21;        %grid fineness
% padsize = 1;          % Pad size is in terms of segment number
% Tseg = 1e-6;          % Thickness of grid segment
% Wseg = 2e-6;          % width of grid segment

T = 25; % (deg C) Temperature

% Package parameters
RPKG = 0.006; % (Ohm)
LPKG = 0.5e-9; % (H)

%% ==================================================
%  ================ BEGIN SIMULATION ================
%  ==================================================
%% Pack objects with inputs
chip.num_gates = Ng;
chip.alpha = alpha;
chip.rent_k = k;
chip.rent_p = p;
chip.num_layers = S;
chip.min_pitch = min_pitch;
chip.gate_pitch = gate_pitch;

chip.area_total = Ach_m2;
chip.chi = chi;
chip.clock_period = Tclk;
chip.logic_activity_factor = a;
chip.Vdd = Vdd;

transistor.gate_length = w_trans;
transistor.oxide_rel_permittivity = eps_ox;
transistor.oxide_thickness = tox;
transistor.leakage_current_per_micron = Ioff;
transistor.capacitance = eps_ox*eps0*w_trans^2/tox;

gate.output_resistance = Ro;
gate.num_transistors = N_trans_per_gate;
gate.capacitance = N_trans_per_gate*transistor.capacitance;

tsv.aspect_ratio = AR_tsv;
tsv.max_area_fraction = Atf_max;
tsv.height = h_tsv_m_thin;

wire.delay_constant = alpha_t;
wire.resistivity = rho_m;
wire.dielectric_epsr = epsr_d;
wire.layers_per_tier = 1;
wire.routing_efficiency = 0.4;
wire.repeater_fraction = [0.5 0.3]; % fraction of optimal repeaters to insert
wire.Beta = [0.25 0.9];
wire.Rc = 0;

simulation.use_joyner = 0;
simulation.redo_wiring_after_repeaters = 0;
simulation.topdown_WLARI = 1;

tic % begin timing
%% TSV number determination
% Inputs:
%   Rent parameters
%   Number of logic gates
%   Number of layers
disp(' ')
disp('Estimating TSV requirements...')
[nt_max nt_tot nt_to nt_through Tacmat] = estimate_tsvs_required(chip.num_gates,chip.num_layers,chip.rent_k,chip.rent_p,chip.alpha);

%% TSV Sizing
% Inputs:
%   Area available
%   Max area for TSVs
%   TSV aspect ratio
disp('Sizing TSVs...')
Ach = Ng/S;
tsv_area_ratio = Atf_max;
[w_tsv_gp h_tsv_gp] = size_tsvs(chip.area_total, tsv.max_area_fraction, nt_max, tsv.aspect_ratio );
%h_tsv_gp = round(h_tsv_gp);
h_tsv_m = h_tsv_gp*gate_pitch;
w_tsv_m = w_tsv_gp*gate_pitch;

tsv.width_m = w_tsv_m;
tsv.height_m = h_tsv_m;
tsv.width_gp = w_tsv_gp;
tsv.height_gp = h_tsv_gp;

%% System determination
% Run WLD + WLA + RI to get power estimate
disp('Generating system...')
%[ iidf l Ln pn pn_orig Cxc Ltot Cn Pdyn Plk Pw Prep Ng_act N_tsvs T_tsvs Atf_act ] = gen_design_old(Ng,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);
[chip power wire repeater tsv] = gen_design(chip,tsv,gate,transistor,wire,simulation);

%Pdens = power.total/chip.area_per_layer_m2;
power.density = power.total/chip.area_total;
psn.pitch_tsv = tsv.pitch_m*100; % [FIX] Power TSVs aren't going to be on the same pitch as signal TSVs

%% System temperature check may go here

%% Power noise estimation
% Inputs:
%   Total power (from previous step)
%   TSV geometry from TSV Sizing step
disp('Evaluating power supply noise...')

psn_iterations = 1;
psn_max = calc_psn(psn,power,chip,tsv,rho_m,mu_m,T);

dispstr = sprintf('\tpsn_runs: %d\tNpads: %d\tpsn_target: %d\tpsn_max: %d',psn_iterations,psn.Npads, psn.noise_target, psn_max);
disp(dispstr)


%% Home in on the best number of power TSVs to use to meet the target
psn.noise = psn_max;
mismatch_norm = psn.noise/psn.noise_target;
psn_iterations_max = 20;
psn_mismatch_tolerance = 0.05; % 5% mismatch is OK
while ( (abs(mismatch_norm-1) > psn_mismatch_tolerance) && (psn_iterations < psn_iterations_max) )
    % more pads -> less noise, and vice versa
    psn.Npads_1d = round(sqrt(psn.Npads*mismatch_norm));
    psn.Npads = psn.Npads_1d^2;
    psn_max = calc_psn(psn,power,chip,tsv,rho_m,mu_m,T);
    psn.noise = psn_max;
    mismatch_norm = psn.noise/psn.noise_target;
    dispstr = sprintf('\tpsn_runs: %d\tNpads: %d\tpsn_target: %d\tpsn_max: %d\tmismatch_norm: %.3g',psn_iterations,psn.Npads, psn.noise_target, psn_max,mismatch_norm);
	disp(dispstr)
end

%% Final report
disp(' ')
disp('Final system parameters:')
repstr = sprintf('\tNg_nom %d \t Ng_act: %d \t Atsv_nom: %.3g \t Atsv_act: %.3g \n\tN_tsvs: %d \t Npads_pow %d \t psn_nom %.4g \t psn_act %.4g', ...
                  Ng, chip.Ng_actual, tsv.max_area_fraction, tsv.actual_area_fraction, tsv.num, psn.Npads, psn.noise_target, psn_max);
disp(repstr)
repstr = sprintf('\th_tsv_um: %.4g \t w_tsv_um: %.4g',h_tsv_m/1e-6,w_tsv_m/1e-6);
disp(repstr);

toc % finish timing

