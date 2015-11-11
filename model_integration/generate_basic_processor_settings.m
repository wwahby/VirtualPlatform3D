function [chip, transistor, gate, tsv, wire, psn, heat] = generate_basic_processor_settings(rent_exp,num_layers,Ng,Ach_mm2,gate_pitch,min_pitch,Vdd,fmax,w_trans)
% simple function to setup standard CPU settings

S = num_layers;

%% constants
eps0 = 8.854e-12; % (F/m) vacuum permittivity


%% Chip parameters
chip.num_gates = Ng;            % (-) Number gates in the system
chip.num_layers = S;            % (-) Number of layers in the 3D stack
chip.area_total = Ach_mm2*1e-6; % (m2) Total chip area
chip.min_pitch = min_pitch;     % (m) Minimum printable pitch (generally contacted gate pitch)
chip.gate_pitch = gate_pitch;   % (m) Actual average gate pitch

chip.fanout = 4;        % average fanout
chip.alpha = chip.fanout/(chip.fanout+1); % input terminal fraction
chip.rent_p = rent_exp;      % rent exponent
chip.rent_k = 3/chip.alpha;  % rent constant
chip.chi = 2/3;         % (-) Conversion factor for point-to-point interconnect length and total net length

chip.clock_period = 1/fmax; % (s) Clock period
chip.logic_activity_factor = 0.1; % (-) Fraction of gates switching at every clock cycle
chip.Vdd = Vdd;        % (V) Supply voltage
chip.temperature = 90;  % (deg C) Temperature guess
chip.thickness_nominal = 50e-6; % (m) Nominal substrate thickness

%% Transistor and gate parameters
transistor.gate_length = w_trans;
transistor.oxide_rel_permittivity = 25; % HfO2
transistor.oxide_thickness = 1.0e-9; %2.5e-9; %1.0e-9 %(nm)
transistor.leakage_current_per_micron = 100e-9; %(A/um)
transistor.leakage_reference_temperature = 390; % (K)
transistor.capacitance_per_micron = 1.25e-15; % (F/um)
transistor.capacitance = transistor.oxide_rel_permittivity*eps0*transistor.gate_length^2/transistor.oxide_thickness;
transistor.subthreshold_swing = .060; % (V/decade at 300K)
transistor.Vt = 0.33; % (V) - Threshold voltage

gate.output_resistance = 8e3;   % (Ohm) Output resistance of a minimum-sized 2in NAND gate
gate.num_transistors = 4;       % (-) number of transistors per average logic gate
gate.capacitance = gate.num_transistors*transistor.capacitance;

%% TSV parameters
tsv.aspect_ratio = 20;          % (-) TSV height / TSV width
tsv.max_area_fraction = 0.10;   % (-) % Maximum fraction of total chip area that TSVs are allowed to consume


%% Wiring parameters
wire.aspect_ratio = 1.8;        % (-) h/w of wires in metal layers
wire.width_fraction = 0.5;     % (-) width/pitch of wires in metal layers

wire.resistivity = 17.2e-9;     % (Ohm*m) Copper wires
wire.barrier_thickness = 0;     % (nm) Thickness of diffusion barrier
wire.barrier_resistivity = 17.2e-9; % (Ohm*m) Resistivity of diffusion barrier
wire.permeability_rel = 1;      % (-) Relative permeability of wiring material
wire.dielectric_epsr = 3.0;     % (-) Relative dielectric constant for wiring ILD -- Low-K dielectric

wire.layers_per_tier = 1;       % (-) Number of metal layers sharing same pitch in each tier
wire.routing_efficiency = [ 0.5 ];  % (-) Fraction of available area that the wire routing tool can actually use
wire.repeater_fraction = [ 0.4 ]; % (-) fraction of optimal repeaters to insert

wire.Beta = [0.9];              % (-v) Fraction of total clock period that a single point-to-point interconnect can consume
wire.Beta_short = 0.25;         % (-) Beta for shortest wiring layers (used for the top down WLARI)

wire.Rc = 0;                    % (-v) Contact resistance between tiers (can be a vector)

wire.use_graphene = 0;          % (-) Allow or disallow graphene use

wire.use_em_resistant_metal = 0;   % (-) Allow or disallow use of electromigration-resistant metals below a specified minimum pitch
wire.min_non_em_width = 50e-9; % (m) If use_em_resistant_metal is set to 1, Cu resistivity will be replaced with wire.alt_resistivity_em below this pitch

wire.alt_resistivity_em = 30e-9;    %(Ohm*m) Resistivity of alternate EM-resistant metal
wire.alt_material_barrier_thickness = 1e-9; % (m) Thickness of diffusion barrier for alternate EM-resistant metal
wire.alt_material_barrier_resistivity = 50e-9; % (Ohm*m) Resistivity of diffusion barrier for alternate EM-resistant metal

wire.repeater_max_area_fraction = 0.2; % (-) Fraction of chip area that can be consumed by repeater/buffer gates
wire.repeater_via_max_area_fraction = 0.05; % (-) Fraction of routable wire area that can be consumed by vias for repeater connections

%% Power supply noise model parameters

psn.noise_fraction = 0.15;          % (-) Maximum power supply noise as a fraction of Vdd
psn.noise_target = 0.1875;          % (V) Acceptable power supply noise
psn.decap_area_fraction = 0.1;      % (Ratio) - Fraction of chip area dedicated to decoupling capacitors
psn.power_tsv_width = 10e-6;        % (m) Width of TSVs used for power delivery. This may be larger than signal TSVs. Setting this to -1 will cause PSN module to use the signal TSV width.

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

% Power TSV determination
psn.mismatch_tolerance = 0.01;      % (-) Allowable normalized deviation from noise target

%% Thermal parameters
%the heat transfer coefficient
% r = 1/(hA); A is the size of top surface area
% the cooling capability of the top heatsink; 20000, 1cm*1cm, means:
% 0.5 W/K
% h = q/dT - q = heat flux (W/m^2)
heat.up = 20000;

% Bottom surface heat transfer coefficient
% This parameter controls the area directly BELOW the bottom chip
% If the interposer is larger than the bottom chip, heat.d controls the
% rest of the area
% Microfluidic heat sinks are assumed to be as large as the chip in the interposer
heat.down = 5;  

% Heat transfer coefficient for the interposer area NOT directly underneath
% the chip(s)
heat.d = 5;

% Side surface heat coefficient, usually near adiabatic
heat.side = 5;

heat.Ta = 298; % ambient temperature

% Alternative settings
heat.r_air = 1/1.825; %K/W for a 1cm^2 HS
heat.r_water = 1/4.63; %K/W for a 1cm^2 HS
heat.A_hs = (1e-2)^2; % 1 cm^2

heat.h_air = 1/(heat.r_air*heat.A_hs);
heat.h_water = 1/(heat.r_water*heat.A_hs);
heat.h_package = 5;

% thermal conductivity of various materials used in thermal module
heat.k_tim = 3;       % TIM
heat.k_chip = 149;      % typically silicon
heat.k_underfill = 0.9; % underfill
heat.k_wires = 400;     % Wires (Typically copper)
heat.k_ild = 1.38;      % ILD
heat.k_microbumps = 60; % Microbumps
heat.k_interposer = 149; % interposer 
heat.k_air = 0.024;      % Air

q_cm2 = 50; % (W/cm2) Top heat sink max heat flux
q = q_cm2*1e4; % (W/m2) Top heat sink max heat flux
dT = 70; % (deg C) Temp difference between chip surface and coolant (air)

heat.up = heat.h_air;
heat.down = heat.h_package;
heat.d = heat.h_package;

heat.interposer_thickness = 200e-6; % (m) Thickness of the interposer below the 3D stack
heat.bump_thickness = 40e-6;        % (m) Microbump thickness (between interposer and bottom chip of 3D stack)
heat.underfill_thickness = 5e-6;    % (m) Thickness of underfill material between each die in the 3D stack
heat.tim_thickness = 5e-6;          % (m) Thickness of thermal interface material between top chip in stack and heat sink
heat.material_IDs = [ 2 9 3];

% If die thickness is thinner than some limit, we're dealing with a
% monolithic 3D stack rather than a conventional 3D stack
heat.monolithic_max_thickness = 9e-6; % (m)
heat.monolithic_intertier_bond_thickness = 0.2e-6;    % (m) Thickness of oxide layer between ILD and next chip
heat.monolithic_material_IDs = [ 2 9 5];

