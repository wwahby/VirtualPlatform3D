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
simulation.force_thickness = 1;



%% 32nm Gulftown (Nehalem arch) single core
Ng_core = 49e6/4; % using transistor number found from other nehalem test case rather than gulftown number, since gulftown number is 50% larger but core is identical
Ach_mm2_core = 19.5;
gate_pitch_core = 2*632e-9; % average gate pitch (sqrt(A_core/Ngates))
min_pitch_core = 112.5e-9; % actual contacted gate pitch
fmax_core = 3.33e9; % 3.33 normal, 3.6 turbo
w_trans = 32e-9;
Vdd_core = 1.30;

%% 32nm Nehalem, SRAM

sram_MB = 6;
Ach_mm2_mem = 32.56;
gate_pitch_mem = 0.804e-6; % average gate pitch (sqrt(A_core/Ngates))
min_pitch_mem = min_pitch_core;
fmax_mem = fmax_core;
Vdd_mem = Vdd_core;

sram_b = sram_MB*2^20;
Ng_mem = sram_b;

%% Queue logic

Ng_queue = 23.2e6/4;
Ach_mm2_queue = 14.53;
gate_pitch_queue = 0.632e-6*2; % average gate pitch (sqrt(A_core/Ngates))
min_pitch_queue = min_pitch_core;
fmax_queue = 733e6;
Vdd_queue = Vdd_core;
w_trans_queue = 32e-9;

%% Memory controller

Ng_con = 55.2e6/4;
Ach_mm2_con = 34.55;
gate_pitch_con = 0.632e-6*2; % average gate pitch (sqrt(A_core/Ngates))
min_pitch_con = min_pitch_core;
fmax_con = fmax_core;
Vdd_con = Vdd_core;
w_trans_con = 32e-9;



%% 
num_layers_per_block = 1;

rent_exp_logic = 0.6;
rent_exp_mem = 0.4;
rent_exp_queue = 0.6;
rent_exp_con = 0.6;

%% define parameters

[core.chip core.transistor core.gate core.tsv core.wire core.psn core.heat] = generate_basic_processor_settings(rent_exp_logic,num_layers_per_block,Ng_core,Ach_mm2_core,gate_pitch_core,min_pitch_core,Vdd_core,fmax_core,w_trans);
[mem.chip mem.transistor mem.gate mem.tsv mem.wire mem.psn mem.heat] = generate_basic_processor_settings(rent_exp_mem,num_layers_per_block,Ng_mem,Ach_mm2_mem,gate_pitch_mem,min_pitch_mem,Vdd_mem,fmax_mem,w_trans);
[queue.chip queue.transistor queue.gate queue.tsv queue.wire queue.psn queue.heat] = generate_basic_processor_settings(rent_exp_queue,num_layers_per_block,Ng_queue,Ach_mm2_queue,gate_pitch_queue,min_pitch_queue,Vdd_queue,fmax_queue,w_trans_queue);
[con.chip con.transistor con.gate con.tsv con.wire con.psn con.heat] = generate_basic_processor_settings(rent_exp_con,num_layers_per_block,Ng_con,Ach_mm2_con,gate_pitch_con,min_pitch_con,Vdd_con,fmax_con,w_trans_queue);


%% Tweak wiring parameters
core.wire.repeater_fraction = [0.4]; % 1 is default from gen_basic_proc_settings
core.wire.routing_efficiency = [0.5]; % 0.4 is default from gen_basic_proc_settings
core.wire.repeater_max_area_fraction = 0.2; % (-) Fraction of chip area that can be consumed by repeater/buffer gates
core.wire.repeater_via_max_area_fraction = 0.05; % (-) Fraction of routable wire area that can be consumed by vias for repeater connections
core.gate.output_resistance = 8e3;
core.transistor.capacitance = 1e-15*1e6*w_trans; % ITRS projection is 1fF/um of gate width. This is an estimate for pMOS transistor capacitance
core.wire.use_graphene = 0;


queue.wire.repeater_fraction = core.wire.repeater_fraction;
queue.wire.routing_efficiency = core.wire.routing_efficiency;
queue.wire.use_graphene = core.wire.use_graphene;
queue.wire.repeater_max_area_fraction = core.wire.repeater_max_area_fraction;
queue.wire.repeater_via_max_area_fraction = core.wire.repeater_via_max_area_fraction;
queue.gate.output_resistane = core.gate.output_resistance;
queue.transistor_capacitance = core.transistor.capacitance;


con.wire.repeater_fraction = core.wire.repeater_fraction;
con.wire.routing_efficiency = core.wire.routing_efficiency;
con.wire.use_graphene = core.wire.use_graphene;
con.wire.repeater_max_area_fraction = core.wire.repeater_max_area_fraction;
con.wire.repeater_via_max_area_fraction = core.wire.repeater_via_max_area_fraction;
con.gate.output_resistane = core.gate.output_resistance;
con.transistor_capacitance = core.transistor.capacitance;

mem.wire.repeater_fraction = core.wire.repeater_fraction;
mem.wire.routing_efficiency = core.wire.routing_efficiency;
mem.wire.use_graphene = core.wire.use_graphene;
mem.wire.repeater_max_area_fraction = core.wire.repeater_max_area_fraction;
mem.wire.repeater_via_max_area_fraction = core.wire.repeater_via_max_area_fraction;
mem.gate.output_resistane = core.gate.output_resistance;
mem.transistor_capacitance = core.transistor.capacitance;


%% calculate block parameters
[core.chip core.power core.tsv core.wire core.repeater core.psn] = codesign_block(core.chip,core.tsv,core.gate,core.transistor,core.wire,core.heat,core.psn,simulation);
[mem.chip mem.power mem.tsv mem.wire mem.repeater mem.psn] = codesign_block(mem.chip,mem.tsv,mem.gate,mem.transistor,mem.wire,mem.heat,mem.psn,simulation);
[queue.chip queue.power queue.tsv queue.wire queue.repeater queue.psn] = codesign_block(queue.chip,queue.tsv,queue.gate,queue.transistor,queue.wire,queue.heat,queue.psn,simulation);
[con.chip con.power con.tsv con.wire con.repeater con.psn] = codesign_block(con.chip,con.tsv,con.gate,con.transistor,con.wire,con.heat,con.psn,simulation);

chip_power_total = 6*core.power.total + mem.power.total + queue.power.total + con.power.total

%% Define block dimensions (m)
core_width = 3.16e-3;
core_height = 6.18e-3;

mem_width = 3*core_width;
mem_height = 3.46e-3;

queue_width = 1.51e-3;
queue_height = 9.64e-3;

con_width = 16.48e-3;
con_height = 2.10e-3;

die1_width = 6*core_width + queue_width;
die1_height = core_height + mem_height + con_height;

cores_x_coords = (0:5) * core_width;
cores_x_coords(4:6) = cores_x_coords(4:6) + queue_width;
cores_y_coords = ones(1,6)*mem_height;

mem1_left_x = 0;
mem1_bot_y = 0;
mem2_left_x = 3*core_width + queue_width;
mem2_bot_y = 0;

queue_x = 3*core_width;
queue_y = 0;

memcon_x = die1_width/2 - con_width/2;
memcon_y = mem_height + core_height;

core_power = core.power.total;
mem_power = mem.power.total;
queue_power = queue.power.total;
con_power = con.power.total;

die1_x_vec = [cores_x_coords mem1_left_x mem2_left_x queue_x memcon_x];
die1_y_vec = [cores_y_coords mem1_bot_y mem2_bot_y queue_y memcon_y];
die1_widths = [ ones(1,6)*core_width mem_width mem_width queue_width con_width];
die1_heights = [ ones(1,6)*core_height mem_height mem_height queue_height con_height];
die1_power_vec = [ ones(1,6)*core_power mem_power mem_power queue_power con_power];

total_chip_width = die1_width;
total_chip_height = die1_height;

map_x_vec = [die1_x_vec]';
map_y_vec = [die1_y_vec]';
map_width = [die1_widths]';
map_height = [die1_heights]';
map_power = [die1_power_vec]';

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
[max_temp temp_vec] = get_stack_temperature(core.chip.num_layers,core.chip.thickness,core.wire,core.tsv,total_chip_width,total_chip_height,package_width,package_height,core.heat,simulation,map,blk_num,power_therm_vec)


%% WLA Validation

intel_pitch = [112.5 112.5 112.5 168.8 225 337.6 450.1 566.5 19400  ]; % Actual intel data
intel_power_total = 130;

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
title('32nm Nehalem (Gulftown)')
