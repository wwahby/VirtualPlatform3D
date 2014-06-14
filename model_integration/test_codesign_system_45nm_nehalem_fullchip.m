% Power and signal codesign
close all
clear all

%% ==================================================
%  ================== BEGIN INPUTS ==================
%  ==================================================

%% Stack parameters
S = 1;

close all
clear all

%% Simulation parameters
simulation.use_joyner = 0;
simulation.redo_wiring_after_repeaters = 0;
simulation.topdown_WLARI = 1; % Use topdown simultaneous WLA and RI (0 = use standard bottom-up optimal WLA, followed by one pass of RI)
simulation.skip_psn_loops = 1; % Skip PSN TSV homing for faster debug
simulation.draw_thermal_map = 0; % Plot thermal profile of each chip
simulation.print_thermal_data = 0; % Output max temp in each layer to console
simulation.separate_wiring_tiers = 0;

%% Logic core parameters
Ach_mm2_core = 32.15;
Ng_core = 33e6/4;
gate_pitch_core = 2*985e-9; % average gate pitch
min_pitch_core = 160e-9; % actual contacted gate pitch
fmax_core = 3.07e9;
w_trans = 45e-9;
Vdd_core = 1.40;

%% Cache parameters
sram_MB = 8;
sram_b = sram_MB*2^20;
Ng_mem = sram_b;
Ach_mm2_mem = 84;
gate_pitch_mem = 1.12e-6; % average gate pitch (sqrt(A_core/Ngates))
min_pitch_mem = 160e-9; % actual contacted gate pitch
fmax_mem = 3.33e9;
Vdd_mem = 1.40;

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

%% 
num_layers_per_block = 1;

rent_exp_logic = 0.6;
rent_exp_mem = 0.4;
rent_exp_gpu = 0.55;

%% define parameters

[core.chip core.transistor core.gate core.tsv core.wire core.psn] = generate_basic_processor_settings(rent_exp_logic,num_layers_per_block,Ng_core,Ach_mm2_core,gate_pitch_core,min_pitch_core,Vdd_core,fmax_core,w_trans);
[mem.chip mem.transistor mem.gate mem.tsv mem.wire mem.psn] = generate_basic_processor_settings(rent_exp_mem,num_layers_per_block,Ng_mem,Ach_mm2_mem,gate_pitch_mem,min_pitch_mem,Vdd_mem,fmax_mem,w_trans);


%% calculate block parameters
[core.chip core.power core.tsv core.wire core.repeater core.psn] = codesign_block(core.chip,core.tsv,core.gate,core.transistor,core.wire,heat,core.psn,simulation);
[mem.chip mem.power mem.tsv mem.wire mem.repeater mem.psn] = codesign_block(mem.chip,mem.tsv,mem.gate,mem.transistor,mem.wire,heat,mem.psn,simulation);

%% Define area
core_width = 4.22e-3;
core_height = 7.62e-3;

mem_width = 18.28e-3;
mem_height = 4.63e-3;

chip_width = 23.29e-3;
chip_height = 12.7e-3;

mem_leftx = 0;
mem_boty = 0;
cores_leftx_vec = [0 core_width (mem_width - 2*core_width) (mem_width - core_width)];
cores_boty_vec = mem_height * [1 1 1 1];
cores_power_vec = core.power.total*ones(1,4);
cores_width_vec = core_width * ones(1,4);
cores_height_vec = core_height * ones(1,4);

map_leftx = [mem_leftx cores_leftx_vec]';
map_boty = [mem_boty cores_boty_vec]';
map_width = [mem_width cores_width_vec]';
map_height = [mem_height cores_height_vec]';
map_power = [mem.power.total cores_power_vec]';

%power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....
map = [map_leftx map_boty map_width map_height map_power];

    %blk_num is for splitting the power maps of each die
    blk_num = [5];
    
%% Thermal module -- Find actual system temperature
chip_power_total = 4*core.power.total + mem.power.total

package_width = 37.5e-3;
package_height = 37.5e-3;
power_therm_vec = chip_power_total;
[max_temp temp_vec] = get_stack_temperature(core.chip.num_layers,core.chip.thickness,core.wire,core.tsv,chip_width,chip_height,package_width,package_height,heat,simulation,map,blk_num,power_therm_vec)
%% WLA Validation

intel_45nm_pitch = [160	160	160	250	280	360	560	810	30500]; % Actual intel data
intel_power_total = 95;
figure(11)
clf
plot(intel_45nm_pitch,'k')
hold on
plot(core.wire.pn*1e9,'r')
xlabel('wiring layer')
ylabel('wire pitch (nm)')
grid on
ylim([ 0 1.1*max([ intel_45nm_pitch(end-1) core.wire.pn(end)*1e9] )]);
fixfigs(11,3,14,12)
title('45nm Nehalem (Lynnfield)')

figure(12)
clf
semilogy(core.wire.wire_area./core.wire.layer_area,'k')
hold on
semilogy(core.wire.via_area_wires./core.wire.layer_area,'b')
semilogy(core.wire.via_area_repeaters./core.wire.layer_area,'r')
fixfigs(12,3,14,12)
