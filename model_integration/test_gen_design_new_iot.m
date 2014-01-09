clear all
close all

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
Ng = 95e6/4;
Ach_mm2 = 11.5;
gate_pitch = 348e-9;
min_pitch = 90e-9;
fmax = 3.0e9;

%% Arbitrarily huge test case
% Ng = 1e9;
% Ach_mm2 = 100;
% gate_pitch = 100e-9;
% min_pitch = 100e-9;
% fmax = 3.0e9;
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
Vdd = 1.0; % (V)

%% constants
eps0 = 8.854e-12; % (F/m) vacuum permittivity

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
%gate.pitch = gate_pitch;

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


%% Corrected distribution
h_tsv_m = h_tsv_m_thin;
tsv.height = h_tsv_m_thin;

[ iidf_old l_old Ln_old pn_old pn_orig_old Cxc_old Ltot_old Cn_old Pdyn_old Plk_old Pw_old Prep_old Ng_act_old N_tsvs_old ] ...
    = gen_design_old(Ng,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,simulation.use_joyner,simulation.redo_wiring_after_repeaters);

chip.num_layers = S;
[chip_new power_new wire_new repeater_new] = gen_design(chip,tsv,gate,transistor,wire,simulation);

%% Plots
figure(1)
clf
plot(repeater_new.num_per_wire)
xlabel('wire length (GP)')
ylabel('Number of repeaters per wire')

figure(2)
clf
plot(repeater_new.size)
xlabel('wire length (GP)')
ylabel('Repeater size (W/Wmin)')

fixfigs(1:2,3,14,12)
