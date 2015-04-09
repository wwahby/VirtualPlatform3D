% Power and signal codesign
% close all
clear all
close all

%% ==================================================
%  ================== BEGIN INPUTS ==================
%  ==================================================

%% Stack parameters
S = 1;


%% 45nm Penryn, entire chip
% Ng = 410e6/4;
% Ach_mm2 = 107;
% gate_pitch = 2*510e-9; % average gate pitch (sqrt(A_core/Ngates))
% min_pitch = 160e-9; % actual contacted gate pitch
% fmax = 3.0e9;
% w_trans = 45e-9;

%% 45nm Penryn, single core
% Ng = 48.6e6/4;
% Ach_mm2 = 21.67;
% gate_pitch = 1.33e-6; % average gate pitch (sqrt(A_core/Ngates))
% min_pitch = 160e-9; % actual contacted gate pitch
% fmax = 3.33e9;
% w_trans = 45e-9;

Ng = 48.6e6/4;
Ach_mm2 = 22.58;
gate_pitch = 1.36e-6; % average gate pitch (sqrt(A_core/Ngates))
min_pitch = 160e-9; % actual contacted gate pitch
fmax = 3.33e9;
w_trans = 45e-9;

%% constants
eps0 = 8.854e-12; % (F/m) vacuum permittivity
%gate_pitch = sqrt( Ach_mm2/1e6/Ng);

rent_exp_logic = 0.6;
num_layers_per_block = 1;
Ng_core = Ng;
Ach_mm2_core = Ach_mm2;
gate_pitch_core = gate_pitch;
min_pitch_core = min_pitch;
Vdd_core = 1.3625;
fmax_core = fmax;

[chip, transistor, gate, tsv, wire, psn, heat] = generate_basic_processor_settings(rent_exp_logic,num_layers_per_block,Ng_core,Ach_mm2_core,gate_pitch_core,min_pitch_core,Vdd_core,fmax_core,w_trans);
%wire.routing_efficiency = [ 0.45 ];  % (-) Fraction of available area that the wire routing tool can actually use
%wire.repeater_fraction = [ 0.4 ]; % (-) fraction of optimal repeaters to insert
wire.routing_efficiency = 0.5;
wire.repeater_fraction = 0.4;

%% Simulation parameters
simulation.use_joyner = 0;
simulation.redo_wiring_after_repeaters = 0;
simulation.topdown_WLARI = 1; % Use topdown simultaneous WLA and RI (0 = use standard bottom-up optimal WLA, followed by one pass of RI)
simulation.skip_psn_loops = 1; % Skip PSN TSV homing for faster debug
simulation.draw_thermal_map = 1; % Plot thermal profile of each chip
simulation.print_thermal_data = 0; % Output max temp in each layer to console
simulation.separate_wiring_tiers = 0;
simulation.force_thickness = 1;



%% Codesign system
tic % begin timing
[core.chip core.power core.tsv core.wire core.repeater core.psn] = codesign_block(chip,tsv,gate,transistor,wire,heat,psn,simulation);
toc % finish timing


%% 45nm Penryn, SRAM
Ng = 302e6/6;
Ach_mm2 = 62;
gate_pitch = 1.11e-6; % average gate pitch (sqrt(A_core/Ngates))
min_pitch = 160e-9; % actual contacted gate pitch
fmax = 3.33e9;
w_trans = 65e-9;

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
chip.Vdd = 1.3625;        % (V) Supply voltage
%chip.temperature = 25;  % (deg C) Temperature
chip.thickness_nominal = 50e-6; % (m) Nominal substrate thickness

tic % begin timing
[mem.chip mem.power mem.tsv mem.wire mem.repeater mem.psn] = codesign_block(chip,tsv,gate,transistor,wire,heat,psn,simulation);
toc % finish timing


%% Define block dimensions (m)
core_width = 6.03e-3;
core_height = 3.74e-3;

mem_width = 6.16e-3;
mem_height = 6.219e-3;

chip_width = 12.36e-3;
chip_height = 8.66e-3;

core1_left_x = 0;
core1_bot_y = (chip_height - 2*core_height)/2;

core2_left_x = 0;
core2_bot_y = core1_bot_y + core_height;

mem_left_x = core_width;
mem_bot_y = (chip_height - mem_height)/2;


%% Set up block map
    %power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....
    map = [mem_left_x     mem_bot_y    mem_width      mem_height      mem.power.total; % L3 Cache
           core1_left_x     core1_bot_y     core_width     core_height     core.power.total;
           core2_left_x     core2_bot_y     core_width     core_height     core.power.total;];
    %blk_num is for splitting the power maps of each die
    blk_num = [3];


simulation.draw_thermal_map = 1; % Plot thermal profile of each chip
simulation.print_thermal_data = 1; % Output max temp in each layer to console

%% Thermal module -- Find actual system temperature
chip_power_total = 2*core.power.total + mem.power.total

package_width = 37.5e-3;
package_height = 37.5e-3;
power_therm_vec = chip_power_total;
[max_temp temp_vec] = get_stack_temperature(core.chip.num_layers,core.chip.thickness,core.wire,core.tsv,chip_width,chip_height,package_width,package_height,heat,simulation,map,blk_num,power_therm_vec)


%% WLA Validation

intel_penryn_pitch = [160	160	160	250	280	360	560	810	30500]; % Actual intel data
intel_power_total = 65;
figure(11)
clf
plot(intel_penryn_pitch,'k')
hold on
plot(core.wire.pn*1e9,'r')
xlabel('wiring layer')
ylabel('wire pitch (nm)')
grid on
ylim([ 0 1.1*intel_penryn_pitch(end-1)]);
fixfigs(11,3,14,12)
title('45nm Core (Penryn/Wolfdale)')

figure(12)
clf
semilogy(core.wire.wire_area./core.wire.layer_area,'k')
hold on
semilogy(core.wire.via_area_wires./core.wire.layer_area,'b')
semilogy(core.wire.via_area_repeaters./core.wire.layer_area,'r')
fixfigs(12,3,14,12)
