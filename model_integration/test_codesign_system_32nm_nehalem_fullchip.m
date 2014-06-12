% Power and signal codesign
% close all
clear all

%% ==================================================
%  ================== BEGIN INPUTS ==================
%  ==================================================

%% Stack parameters
S = 1;

%% 45nm Nehalem, single core
%Ng = 48e6/4;
Ng = 39e6/4;
Ach_mm2 = 18.3;
gate_pitch = 2*693e-9; % average gate pitch (sqrt(A_core/Ngates))
min_pitch = 112.5e-9; % actual contacted gate pitch
fmax = 3.6e9;
w_trans = 32e-9;


%% constants
eps0 = 8.854e-12; % (F/m) vacuum permittivity
%gate_pitch = sqrt( Ach_mm2/1e6/Ng);

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
chip.Vdd = 1.40;        % (V) Supply voltage
chip.temperature = 70;  % (deg C) Temperature
chip.thickness_nominal = 50e-6; % (m) Nominal substrate thickness

%% Transistor and gate parameters
transistor.gate_length = w_trans;
transistor.oxide_rel_permittivity = 25; % HfO2
transistor.oxide_thickness = 1e-9;
transistor.leakage_current_per_micron = 10e-9; %(A/um)
transistor.capacitance = transistor.oxide_rel_permittivity*eps0*transistor.gate_length^2/transistor.oxide_thickness;
transistor.subthreshold_swing = .060; % (V/decade at 300K)
transistor.Vt = 0.25; % (V) - Threhsold voltage

gate.output_resistance = 5e3;   % (Ohm) Output resistance of a minimum-sized 2in NAND gate
gate.num_transistors = 4;       % (-) number of transistors per average logic gate
gate.capacitance = gate.num_transistors*transistor.capacitance;

%% TSV parameters
tsv.aspect_ratio = 20;          % (-) TSV height / TSV width
tsv.max_area_fraction = 0.10;   % (-) % Maximum fraction of total chip area that TSVs are allowed to consume


%% Wiring parameters
wire.aspect_ratio = 1.8;        % (-) h/w of wires in metal layers
wire.width_fraction = 0.5;     % (-) width/pitch of wires in metal layers
wire.resistivity = 17.2e-9;     % (Ohm*m) Copper wires
wire.permeability_rel = 1;      % (-) Relative permeability of wiring material
wire.dielectric_epsr = 3.0;     % (-) Relative dielectric constant for wiring ILD -- Low-K dielectric
wire.layers_per_tier = 1;       % (-) Number of metal layers sharing same pitch in each tier
wire.routing_efficiency = [0.38 ];  % (-) Fraction of available area that the wire routing tool can actually use
wire.repeater_fraction = [ 0.5 ]; % (-) fraction of optimal repeaters to insert
wire.Beta = [0.9];              % (-v) Fraction of total clock period that a single point-to-point interconnect can consume
wire.Beta_short = 0.25;         % (-) Beta for shortest wiring layers (used for the top down WLARI)
wire.Rc = 0;                    % (-v) Contact resistance between tiers (can be a vector)
wire.use_graphene = 0;

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

% Power TSV determination
psn.mismatch_tolerance = 0.05;      % (-) Allowable normalized deviation from noise target

mu_m = 1.257e-6;      %copper permeability

% decap = 0.1; % (Ratio) - Fraction of chip area dedicated to decoupling capacitors
% Npad1d = 32;
% Npad = Npad1d^2;         %total number of power or ground pads
% Ngrid = 21*21;        %grid fineness
% padsize = 1;          % Pad size is in terms of segment number
% Tseg = 1e-6;          % Thickness of grid segment
% Wseg = 2e-6;          % width of grid segment

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
q_cm2 = 50; % (W/cm2) Top heat sink max heat flux
q = q_cm2*1e4; % (W/m2) Top heat sink max heat flux
dT = 70; % (deg C) Temp difference between chip surface and coolant (air)
heat.up = q/dT;
heat.down = 2*heat.up;
heat.d = heat.down;

%% Simulation parameters
simulation.use_joyner = 0;
simulation.redo_wiring_after_repeaters = 0;
simulation.topdown_WLARI = 1; % Use topdown simultaneous WLA and RI (0 = use standard bottom-up optimal WLA, followed by one pass of RI)
simulation.skip_psn_loops = 1; % Skip PSN TSV homing for faster debug
simulation.draw_thermal_map = 0; % Plot thermal profile of each chip
simulation.print_thermal_data = 0; % Output max temp in each layer to console
simulation.separate_wiring_tiers = 1;



%% Codesign system
tic % begin timing
[core.chip core.power core.tsv core.wire core.repeater core.psn] = codesign_system(chip,tsv,gate,transistor,wire,heat,psn,simulation);
toc % finish timing


%% 32nm Nehalem, SRAM

sram_MB = 4;
Ach_mm2 = 20.8;
gate_pitch = 0.788e-6; % average gate pitch (sqrt(A_core/Ngates))

sram_b = sram_MB*2^20;
Ng = sram_b;


chip.num_gates = Ng;            % (-) Number gates in the system
chip.num_layers = S;            % (-) Number of layers in the 3D stack
chip.area_total = Ach_mm2*1e-6; % (m2) Total chip area
chip.min_pitch = min_pitch;     % (m) Minimum printable pitch (generally contacted gate pitch)
chip.gate_pitch = gate_pitch;   % (m) Actual average gate pitch

chip.fanout = 4;        % average fanout
chip.alpha = chip.fanout/(chip.fanout+1); % input terminal fraction
chip.rent_p = 0.4;      % rent exponent
chip.rent_k = 3/chip.alpha;  % rent constant
chip.chi = 2/3;         % (-) Conversion factor for point-to-point interconnect length and total net length

chip.clock_period = 1/fmax; % (s) Clock period
chip.logic_activity_factor = 0.01; % (-) Fraction of gates switching at every clock cycle
chip.Vdd = 1.40;        % (V) Supply voltage
%chip.temperature = 25;  % (deg C) Temperature
chip.thickness_nominal = 50e-6; % (m) Nominal substrate thickness

tic % begin timing
[mem.chip mem.power mem.tsv mem.wire mem.repeater mem.psn] = codesign_system(chip,tsv,gate,transistor,wire,heat,psn,simulation);
toc % finish timing


%% 32nm Nehalem -- 45nm GPU in MCP

Ng = 114e6/4;
Ach_mm2 = 114;
gate_pitch = 0.802e-6*2; % average gate pitch (sqrt(A_core/Ngates))

chip.num_gates = Ng;            % (-) Number gates in the system
chip.num_layers = S;            % (-) Number of layers in the 3D stack
chip.area_total = Ach_mm2*1e-6; % (m2) Total chip area
chip.min_pitch = min_pitch;     % (m) Minimum printable pitch (generally contacted gate pitch)
chip.gate_pitch = gate_pitch;   % (m) Actual average gate pitch

chip.fanout = 4;        % average fanout
chip.alpha = chip.fanout/(chip.fanout+1); % input terminal fraction
chip.rent_p = 0.4;      % rent exponent
chip.rent_k = 3/chip.alpha;  % rent constant
chip.chi = 2/3;         % (-) Conversion factor for point-to-point interconnect length and total net length

chip.clock_period = 1/fmax; % (s) Clock period
chip.logic_activity_factor = 0.1; % (-) Fraction of gates switching at every clock cycle
chip.Vdd = 1.40;        % (V) Supply voltage
%chip.temperature = 25;  % (deg C) Temperature
chip.thickness_nominal = 50e-6; % (m) Nominal substrate thickness

tic % begin timing
[gpu.chip gpu.power gpu.tsv gpu.wire gpu.repeater gpu.psn] = codesign_system(chip,tsv,gate,transistor,wire,heat,psn,simulation);


%%
num_cores = 2;
chip_power_total = num_cores*core.power.total + mem.power.total + gpu.power.total

%% WLA Validation

intel_pitch = [112.5 112.5 112.5 168.8 225 337.6 450.1 566.5 19400  ]; % Actual intel data
intel_power_total = 73;

figure(1)
clf
plot(intel_pitch,'k')
hold on
plot(core.wire.pn*1e9,'r')
xlabel('wiring layer')
ylabel('wire pitch (nm)')
grid on
ylim([ 0 1.1*max([ intel_pitch(end-1) core.wire.pn(end)*1e9] )]);
fixfigs(1,3,14,12)
title('32nm Nehalem (Clarkdale)')

figure(2)
clf
semilogy(core.wire.wire_area./core.wire.layer_area,'k')
hold on
semilogy(core.wire.via_area_wires./core.wire.layer_area,'b')
semilogy(core.wire.via_area_repeaters./core.wire.layer_area,'r')
fixfigs(2,3,14,12)

figure(3)
clf
bar([chip_power_total intel_power_total])