% Power and signal codesign
% close all
clear all

%% ==================================================
%  ================== BEGIN INPUTS ==================
%  ==================================================

%% Stack parameters
S = 1;

%% Simulation parameters
simulation.use_joyner = 0;
simulation.redo_wiring_after_repeaters = 0;
simulation.topdown_WLARI = 1; % Use topdown simultaneous WLA and RI (0 = use standard bottom-up optimal WLA, followed by one pass of RI)
simulation.skip_psn_loops = 1; % Skip PSN TSV homing for faster debug
simulation.draw_thermal_map = 1; % Plot thermal profile of each chip
simulation.print_thermal_data = 0; % Output max temp in each layer to console
simulation.separate_wiring_tiers = 1;



%% 32nm Nehalem, single core
%Ng = 48e6/4;
Ng_core = 39e6/4;
Ach_mm2_core = 18.3;
gate_pitch_core = 2*693e-9; % average gate pitch (sqrt(A_core/Ngates))
min_pitch_core = 112.5e-9; % actual contacted gate pitch
fmax_core = 3.6e9;
w_trans = 32e-9;
Vdd_core = 1.40;


%% 32nm Nehalem, SRAM

sram_MB = 4;
Ach_mm2_mem = 20.8;
gate_pitch_mem = 0.788e-6; % average gate pitch (sqrt(A_core/Ngates))
min_pitch_mem = min_pitch_core;
fmax_mem = fmax_core;
Vdd_mem = Vdd_core;

sram_b = sram_MB*2^20;
Ng_mem = sram_b;

%% 32nm Nehalem -- 45nm GPU in MCP

Ng_gpu = 114e6/4;
Ach_mm2_gpu = 73; %114;
gate_pitch_gpu = 0.802e-6*2; % average gate pitch (sqrt(A_core/Ngates))
min_pitch_gpu = 160e-9;
fmax_gpu = 733e6;
Vdd_gpu = Vdd_core;
w_trans_gpu = 45e-9;

%% 32nm Nehalem -- 45nm PCIE Controller in MCP

Ng_con = 47.1e6/4;
Ach_mm2_con = 30.3;
gate_pitch_con = 0.802e-6*2; % average gate pitch (sqrt(A_core/Ngates))
min_pitch_con = 160e-9;
fmax_con = 900e6;
Vdd_con = Vdd_core;
w_trans_con = 45e-9;



%% 
num_layers_per_block = 1;

rent_exp_logic = 0.6;
rent_exp_mem = 0.4;
rent_exp_gpu = 0.6;
rent_exp_con = 0.6;

%% define parameters

[core.chip core.transistor core.gate core.tsv core.wire core.psn] = generate_basic_processor_settings(rent_exp_logic,num_layers_per_block,Ng_core,Ach_mm2_core,gate_pitch_core,min_pitch_core,Vdd_core,fmax_core,w_trans);
[mem.chip mem.transistor mem.gate mem.tsv mem.wire mem.psn] = generate_basic_processor_settings(rent_exp_mem,num_layers_per_block,Ng_mem,Ach_mm2_mem,gate_pitch_mem,min_pitch_mem,Vdd_mem,fmax_mem,w_trans);
[gpu.chip gpu.transistor gpu.gate gpu.tsv gpu.wire gpu.psn] = generate_basic_processor_settings(rent_exp_gpu,num_layers_per_block,Ng_gpu,Ach_mm2_gpu,gate_pitch_gpu,min_pitch_gpu,Vdd_gpu,fmax_gpu,w_trans_gpu);
[con.chip con.transistor con.gate con.tsv con.wire con.psn] = generate_basic_processor_settings(rent_exp_con,num_layers_per_block,Ng_con,Ach_mm2_con,gate_pitch_con,min_pitch_con,Vdd_con,fmax_con,w_trans_gpu);


%% Tweak wiring parameters
core.wire.repeater_fraction = [0.3]; % 1 is default from gen_basic_proc_settings
core.wire.routing_efficiency = [0.4]; % 0.4 is default from gen_basic_proc_settings
core.gate.output_resistance = 8e3; % Ohm
core.wire.use_graphene = 0;

gpu.wire.repeater_fraction = core.wire.repeater_fraction;
gpu.wire.routing_efficiency = core.wire.routing_efficiency;
gpu.wire.use_graphene = core.wire.use_graphene;
gpu.gate.output_resistance = 8e3;

con.wire.repeater_fraction = core.wire.repeater_fraction;
con.wire.routing_efficiency = core.wire.routing_efficiency;
con.wire.use_graphene = core.wire.use_graphene;
con.gate.output_resistance = 8e3;

mem.wire.repeater_fraction = core.wire.repeater_fraction;
mem.wire.routing_efficiency = core.wire.routing_efficiency;
mem.wire.use_graphene = core.wire.use_graphene;
mem.gate.output_resistance = 8e3;


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
% q_cm2 = 50; % (W/cm2) Top heat sink max heat flux
% q = q_cm2*1e4; % (W/m2) Top heat sink max heat flux
% dT = 70; % (deg C) Temp difference between chip surface and coolant (air)
% heat.up = q/dT;
% heat.down = 5;
% heat.d = 5;
% heat.side = 5;


%% calculate block parameters
[core.chip core.power core.tsv core.wire core.repeater core.psn] = codesign_block(core.chip,core.tsv,core.gate,core.transistor,core.wire,heat,core.psn,simulation);
[mem.chip mem.power mem.tsv mem.wire mem.repeater mem.psn] = codesign_block(mem.chip,mem.tsv,mem.gate,mem.transistor,mem.wire,heat,mem.psn,simulation);
[gpu.chip gpu.power gpu.tsv gpu.wire gpu.repeater gpu.psn] = codesign_block(gpu.chip,gpu.tsv,gpu.gate,gpu.transistor,gpu.wire,heat,gpu.psn,simulation);
[con.chip con.power con.tsv con.wire con.repeater con.psn] = codesign_block(con.chip,con.tsv,con.gate,con.transistor,con.wire,heat,con.psn,simulation);

chip_power_total = 2*core.power.total + mem.power.total + gpu.power.total + con.power.total

%% Define block dimensions (m)
core_width = 3.07e-3;
core_height = 5.96e-3;

mem_width = 6.02e-3;
mem_height = 3.46e-3;

die1_width = 8.49e-3;
die1_height = 9.55e-3;

die2_width = 11.45e-3;
die2_height = 15.46e-3;

con_width = 5.71e-3;
con_height = 5.31e-3;

gpu_widths = [0.94 8.15 0.88 3.11 1.54]*1e-3;
gpu_heights = [1.78 6.43 1.33 5.40 0.82]*1e-3;

gpu_areas = gpu_heights .* gpu_widths;
gpu_area = sum(gpu_areas);
gpu_rel_areas = gpu_areas/gpu_area;
gpu_rel_powers = gpu.power.total * gpu_rel_areas;

gpu1_leftx = die2_width - sum(gpu_widths(1:3));
gpu2_leftx = gpu1_leftx + gpu_widths(1);
gpu3_leftx = gpu2_leftx + gpu_widths(2);
gpu4_leftx = die2_width - gpu_widths(3) - gpu_widths(4);
gpu5_leftx = die2_width - gpu_widths(3) - gpu_widths(5);

gpu1_boty = gpu_heights(2) - gpu_heights(1);
gpu2_boty = 0;
gpu3_boty = 0;
gpu4_boty = gpu_heights(2);
gpu5_boty = gpu4_boty + gpu_heights(4);

gpu_left_x_coords = [gpu1_leftx gpu2_leftx gpu3_leftx gpu4_leftx gpu5_leftx];
gpu_bot_y_coords = [gpu1_boty gpu2_boty gpu3_boty gpu4_boty gpu5_boty];

mem_io_height = 1.8e-3;
gpu_bot_y_coords = gpu_bot_y_coords + mem_io_height; % stand off from bottom of chip
gpu_power_vec = gpu_rel_powers;

con_leftx = gpu_left_x_coords(1);
con_boty = gpu_bot_y_coords(1) + gpu_heights(1);
con_power = con.power.total;

die2_leftx = [con_leftx gpu_left_x_coords];
die2_boty = [con_boty gpu_bot_y_coords];
die2_widths = [con_width gpu_widths];
die2_heights = [con_height gpu_heights];
die2_power = [con_power gpu_power_vec];


cores_x_coords = [0 core_width];
cores_y_coords = [0 0];

mem_left_x = 0;
mem_bot_y = (die1_height - mem_height);

die1_x_vec = [cores_x_coords mem_left_x];
die1_y_vec = [cores_y_coords mem_bot_y];

% shift cores over
standoff_width = 1e-3;
die1_x_vec = die1_x_vec + die2_width + standoff_width; % 1mm standoff between chips
die1_y_vec = die1_y_vec + (die2_height - die1_height)/2;

die1_widths = [core_width core_width mem_width];
die1_heights = [core_height core_height mem_height];

die1_power_vec = [core.power.total core.power.total mem.power.total];

total_chip_width = die1_width + die2_width + standoff_width;
total_chip_height = max(die1_width,die2_width);



map_x_vec = [die2_leftx die1_x_vec]';
map_y_vec = [die2_boty die1_y_vec]';
map_width = [die2_widths die1_widths]';
map_height = [die2_heights die1_heights]';
map_power = [die2_power die1_power_vec]';

%% Set up block map
    %power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....

    map = [map_x_vec map_y_vec map_width map_height map_power];
    %blk_num is for splitting the power maps of each die
    blk_num = [length(map_x_vec)];


%% Thermal module -- Find actual system temperature

package_width = 37.5e-3;
package_height = 37.5e-3;
power_therm_vec = chip_power_total;
[max_temp temp_vec] = get_stack_temperature(core.chip.num_layers,core.chip.thickness,core.wire,core.tsv,total_chip_width,total_chip_height,package_width,package_height,heat,simulation,map,blk_num,power_therm_vec)


%% WLA Validation

intel_pitch = [112.5 112.5 112.5 168.8 225 337.6 450.1 566.5 19400  ]; % Actual intel data
intel_power_total = 73;

figure(11)
clf
plot(intel_pitch,'k')
hold on
plot(core.wire.pn*1e9,'r')
xlabel('wiring layer')
ylabel('wire pitch (nm)')
grid on
ylim([ 0 1.1*max([ intel_pitch(end-1) core.wire.pn(end)*1e9] )]);
fixfigs(11,3,14,12)
title('32nm Nehalem (Clarkdale)')

figure(12)
clf
semilogy(core.wire.wire_area./core.wire.layer_area,'k')
hold on
semilogy(core.wire.via_area_wires./core.wire.layer_area,'b')
semilogy(core.wire.via_area_repeaters./core.wire.layer_area,'r')
fixfigs(12,3,14,12)