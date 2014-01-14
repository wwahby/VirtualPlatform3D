% Power and signal codesign
close all
clear all

%% ==================================================
%  ================== BEGIN INPUTS ==================
%  ==================================================

%% Stack parameters
S = 2;

%% 65nm Merom, entire chip
% Ng = 291e6/4;
% Ach_mm2 = 143;
% gate_pitch = 700e-9; % average gate pitch (sqrt(A_core/Ngates))
% min_pitch = 220e-9; % actual contacted gate pitch
% fmax = 3.0e9;
% w_trans = 65e-9;

%% 65nm Merom, single core, wrong area
% Ng = 10e6;
% Ach_mm2 = 16;
% gate_pitch = 630e-9; % average gate pitch (sqrt(A_core/Ngates))
% min_pitch = 220e-9; % actual contacted gate pitch
% fmax = 3.0e9;
% w_trans = 65e-9;

%% 65nm Merom, single core
% Ng = 10e6;
% Ach_mm2 = 27;
% gate_pitch = 820e-9; % average gate pitch (sqrt(A_core/Ngates))
% min_pitch = 220e-9; % actual contacted gate pitch
% fmax = 3.0e9;
% w_trans = 65e-9;

%% 32nm Sandy Bridge, entire chip
% Ng = 2.27e9/4;
% Ach_mm2 = 435;
% gate_pitch = 875e-9;
% min_pitch = 112.5e-9;
% fmax = 3.6e9;
% w_trans = 32e-9;

%% 32nm Sandy Bridge, one core
% Ng = 107e6/4;
% Ach_mm2 = 20;
% gate_pitch = 435e-9;
% min_pitch = 112.5e-9;
% fmax = 3.6e9;
% w_trans = 32e-9;

%% 22nm Ivy Bridge, entire chip
% Ng = 2.86e9/4;
% Ach_mm2 = 346.5;
% gate_pitch = 348e-9;
% min_pitch = 90e-9;
% fmax = 3.0e9;
% w_trans = 80e-9;

%% 22nm Ivy Bridge, one core
Ng = 95e6/4;
Ach_mm2 = 11.5;
gate_pitch = 348e-9;
min_pitch = 90e-9;
fmax = 3.0e9;
w_trans = 80e-9;

%% Arbitrarily huge test case
% Ng = 1e9;
% Ach_mm2 = 100;
% gate_pitch = 100e-9;
% min_pitch = 100e-9;
% fmax = 3.0e9;
% w_trans = 25e-9;
% %% Chip descriptors
% 
% Ach_m2 = Ach_mm2*1e-6;
% 
% % gate parameters
% eps_ox = 25; % HfO2
% tox = 1e-9;
% 
% N_trans_per_gate = 4;
% Ioff = 10e-9; %(A/um)
% 
% % Tsv parameters
% Atf_max = 0.10; % maximum allowable TSV area, as a fraction of total chip area
% h_tsv_m_thin = 10e-6;
% h_tsv_m_thick = 300e-6;
% AR_tsv = 20;
% 
% % Rent parameters
% p = 0.6; % rent exponent
% fo = 4; % avg fanout
% alpha = fo/(fo+1); % input terminal fraction
% k = 3/alpha; %rent constant
% 
% % Wiring parameters
% chi = 2/3;
% rho_m = 17.2e-9; % Cu
% epsr_d = 3.0; % Low-k dielectric
% Tclk = 1/fmax; % (s)
% alpha_t = 1.1*6.2;
% 
% % Repeater parameters
% Ro = 1e3; % (Ohm) Gate output resistance
% 
% % Power parameters
% a = 0.1; % logic activity factor
% Vdd = 1.25; % (V)

%% constants
eps0 = 8.854e-12; % (F/m) vacuum permittivity

%% ==================================================
%  ================ BEGIN SIMULATION ================
%  ==================================================
%% Chip parameters
chip.num_gates = Ng;            % (-) Number gates in the system
chip.num_layers = S;            % (-) Number of layers in the 3D stack
chip.area_total = Ach_mm2*1e-6; % (m2) Total chip area
chip.min_pitch = min_pitch;     % (m) Minimum printable pitch (generally contacted gate pitch)
chip.gate_pitch = gate_pitch;   % (m) Actual average gate pitch

chip.fanout = 4;        % average fanout
chip.alpha = chip.fanout/(chip.fanout+1); % input terminal fraction
chip.rent_p = 0.6;      % rent exponent
chip.rent_k = 3/chip.alpha;  % rent constant
chip.chi = 2/3;         % (-) Conversion factor for point-to-point interconnect length and total net length

chip.clock_period = 1/fmax; % (s) Clock period
chip.logic_activity_factor = 0.1; % (-) Fraction of gates switching at every clock cycle
chip.Vdd = 1.25;        % (V) Supply voltage
chip.temperature = 25;  % (deg C) Temperature

%% Transistor and gate parameters
transistor.gate_length = w_trans;
transistor.oxide_rel_permittivity = 25; % HfO2
transistor.oxide_thickness = 1e-9;
transistor.leakage_current_per_micron = 10e-9; %(A/um)
transistor.capacitance = transistor.oxide_rel_permittivity*eps0*transistor.gate_length^2/transistor.oxide_thickness;

gate.output_resistance = 1e3;   % (Ohm) Output resistance of a minimum-sized 2in NAND gate
gate.num_transistors = 4;       % (-) number of transistors per average logic gate
gate.capacitance = gate.num_transistors*transistor.capacitance;

%% TSV and wire parameters
tsv.aspect_ratio = 20;          % (-) TSV height / TSV width
tsv.max_area_fraction = 0.10;   % (-) % Maximum fraction of total chip area that TSVs are allowed to consume
tsv.height = 50e-6;             % (m) TSV height

wire.delay_constant = 1.1*6.2;  % (-) Multiplicative constant determining average wire delay (1.1 is to get more than 50% delay, 6.2 comes from average wire capacitance -- see Venkatesan)
wire.resistivity = 17.2e-9;     % (Ohm*m) Copper wires
wire.permeability_rel = 1;      % (-) Relative permeability of wiring material
wire.dielectric_epsr = 3.0;     % (-) Relative dielectric constant for wiring ILD -- Low-K dielectric
wire.layers_per_tier = 1;       % (-) Number of metal layers sharing same pitch in each tier
wire.routing_efficiency = 0.4;  % (-) Fraction of available area that the wire routing tool can actually use
wire.repeater_fraction = [0.5 0.3]; % (-) fraction of optimal repeaters to insert
wire.Beta = [0.9];              % (-v) Fraction of total clock period that a single point-to-point interconnect can consume
wire.Rc = 0;                    % (-v) Contact resistance between tiers (can be a vector)

%% Power supply noise model parameters

psn.noise_fraction = 0.15;          % (-) Maximum power supply noise as a fraction of Vdd
psn.noise_target = 0.1875;          % (V) Acceptable power supply noise
psn.decap_area_fraction = 0.1;      % (Ratio) - Fraction of chip area dedicated to decoupling capacitors

% Power connection parameters
psn.Npads_1d = 100;                  % (-) Number of pads from one side of chip to the other
psn.Npads = psn.Npads_1d^2;         % (-) total number of power pads (does not include ground pads)
psn.Ngrid = 21*21;                  % (-) grid fineness
psn.pad_size = 1;                   % (-) Pad size is in terms of segment number
psn.segment_thickness = 1e-6;       % (m) Thickness of a power grid segment
psn.segment_width = 2e-6;           % (m) width of a power grid segment

% Package parameters
psn.package_resistance = 0.006;     % (Ohm) Resistance per pad on the package
psn.package_inductance = 0.5e-9;    % (H) Inductance per package pad

mu_m = 1.257e-6;      %copper permeability

% decap = 0.1; % (Ratio) - Fraction of chip area dedicated to decoupling capacitors
% Npad1d = 32;
% Npad = Npad1d^2;         %total number of power or ground pads
% Ngrid = 21*21;        %grid fineness
% padsize = 1;          % Pad size is in terms of segment number
% Tseg = 1e-6;          % Thickness of grid segment
% Wseg = 2e-6;          % width of grid segment

%% Simulation parameters
simulation.use_joyner = 0;
simulation.redo_wiring_after_repeaters = 0;
simulation.topdown_WLARI = 1; % Use topdown simultaneous WLA and RI (0 = use standard bottom-up optimal WLA, followed by one pass of RI)

%% Codesign system
tic % begin timing
[chip power tsv wire repeater psn] = codesign_system(chip,tsv,gate,transistor,wire,psn,simulation);
toc % finish timing

