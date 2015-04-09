close all

%% Simulation parameters
simulation.use_joyner = 0;
simulation.redo_wiring_after_repeaters = 0;
simulation.topdown_WLARI = 1; % Use topdown simultaneous WLA and RI (0 = use standard bottom-up optimal WLA, followed by one pass of RI)
simulation.skip_psn_loops = 1; % Skip PSN TSV homing for faster debug
simulation.draw_thermal_map = 0; % Plot thermal profile of each chip
simulation.print_thermal_data = 0; % Output max temp in each layer to console
simulation.separate_wiring_tiers = 1; % 1 = calculate wire pitch for wiring tiers between EACH logic plane
simulation.force_thickness = 1;

%% Logic core parameters

% Ng_core = 86e6/4;
% Ach_mm2_core = 18.5;
% gate_pitch_core = 932e-9;
% min_pitch_core = 112.5e-9;
% fmax_core = 3.5e9; % 3.5GHz normal, 3.9GHz turbo
% w_trans = 32e-9;
% Vdd_core = 1.25;

Ng_core = 86e6/4;
Ach_mm2_core = 18.5;
gate_pitch_core = 465e-9*2;
min_pitch_core = 112.5e-9;
fmax_core = 3.5e9;
w_trans = 32e-9;
Vdd_core = 1.25;

%% Cache parameters
sram_MB = 8;
sram_b = sram_MB*2^20;
Ng_mem = sram_b;
Ach_mm2_mem = 46;
gate_pitch_mem = 0.828e-6; % average gate pitch (sqrt(A_core/Ngates))
min_pitch_mem = 112.5e-9; % actual contacted gate pitch
fmax_mem = 3.5e9;
Vdd_mem = 1.25;

%% GPU parameters
Ng_gpu = 2.05e8/4;
Ach_mm2_gpu = 44.5;
gate_pitch_gpu = 0.931e-6; % average gate pitch (sqrt(A_core/Ngates))
min_pitch_gpu = 112.5e-9; % actual contacted gate pitch
fmax_gpu = 1.35e9; % max clock
%fmax_gpu = 850e6; % base clock 
Vdd_gpu = 1.25;
%% 
num_layers_per_block = 1;

rent_exp_logic = 0.6;
rent_exp_mem = 0.4;
rent_exp_gpu = 0.50;

%% define parameters

[core.chip core.transistor core.gate core.tsv core.wire core.psn core.heat] = generate_basic_processor_settings(rent_exp_logic,num_layers_per_block,Ng_core,Ach_mm2_core,gate_pitch_core,min_pitch_core,Vdd_core,fmax_core,w_trans);
[mem.chip mem.transistor mem.gate mem.tsv mem.wire mem.psn mem.heat] = generate_basic_processor_settings(rent_exp_mem,num_layers_per_block,Ng_mem,Ach_mm2_mem,gate_pitch_mem,min_pitch_mem,Vdd_mem,fmax_mem,w_trans);
[gpu.chip gpu.transistor gpu.gate gpu.tsv gpu.wire gpu.psn gpu.heat] = generate_basic_processor_settings(rent_exp_gpu,num_layers_per_block,Ng_gpu,Ach_mm2_gpu,gate_pitch_gpu,min_pitch_gpu,Vdd_gpu,fmax_gpu,w_trans);

%% Tweak wiring parameters
core.wire.repeater_fraction = [0.4]; % 1 is default from gen_basic_proc_settings
core.wire.routing_efficiency = [0.5]; % 0.4 is default from gen_basic_proc_settings
core.wire.repeater_max_area_fraction = 0.2;
core.wire.repeater_via_max_area_fraction = 0.05;
core.gate.output_resistance = 8e3; % Ohm
core.transistor.capacitance = 1e-15*1e6*3*w_trans; % ITRS projection is 1fF/um of gate width. This is an estimate for pMOS transistor capacitance


gpu.wire.repeater_fraction = core.wire.repeater_fraction;
gpu.wire.routing_efficiency = core.wire.routing_efficiency;
gpu.gate.output_resistance = core.gate.output_resistance; % Ohm
gpu.wire.repeater_max_area_fraction = core.wire.repeater_max_area_fraction;
gpu.wire.repeater_via_max_area_fraction = core.wire.repeater_via_max_area_fraction;


mem.wire.repeater_fraction = core.wire.repeater_fraction;
mem.wire.routing_efficiency = core.wire.routing_efficiency;
mem.gate.output_resistance = core.gate.output_resistance; % Ohm
mem.wire.repeater_max_area_fraction = core.wire.repeater_max_area_fraction;
mem.wire.repeater_via_max_area_fraction = core.wire.repeater_via_max_area_fraction;

%% calculate block parameters
[core.chip core.power core.tsv core.wire core.repeater core.psn] = codesign_block(core.chip,core.tsv,core.gate,core.transistor,core.wire,core.heat,core.psn,simulation);
[mem.chip mem.power mem.tsv mem.wire mem.repeater mem.psn] = codesign_block(mem.chip,mem.tsv,mem.gate,mem.transistor,mem.wire,mem.heat,mem.psn,simulation);
[gpu.chip gpu.power gpu.tsv gpu.wire gpu.repeater gpu.psn] = codesign_block(gpu.chip,gpu.tsv,gpu.gate,gpu.transistor,gpu.wire,gpu.heat,gpu.psn,simulation);


%% Define block dimensions (m)
core_width = 3.302e-3;
core_height = 5.65e-3;

mem_width = 13.24e-3;
mem_height = 3.48e-3;

gpu_width = 4.89e-3;
gpu_height = 9.12e-3;

chip_width = 21.08e-3;
chip_height = 10.25e-3;

mem_io_width = 13.49e-3;
mem_io_height = 1.16e-3;

sys_agent_width = 3.01e-3;
sys_agent_height = 9.11e-3;

core_bot_y = mem_io_height + mem_height;
core1_left_x = gpu_width;
core2_left_x = core1_left_x + core_width;
core3_left_x = core2_left_x + core_width;
core4_left_x = core3_left_x + core_width;

gpu_left_x = 0;
gpu_bot_y = mem_io_height;

mem_left_x = gpu_width;
mem_bot_y = mem_io_height;

mem_io_left_x = chip_width - mem_io_width;
mem_io_bot_y = 0;

sys_agent_left_x = gpu_width + mem_width;
sys_agent_bot_y = mem_io_height;


%% Set up block map
    %power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....
    map = [mem_io_left_x   mem_io_bot_y    mem_io_width   mem_io_height     0; % memory IO
           sys_agent_left_x   sys_agent_bot_y    sys_agent_width   sys_agent_height     0; % system agent + memory controller
           gpu_left_x      gpu_bot_y       gpu_width      gpu_height      gpu.power.total; % GPU
           mem_left_x     mem_io_height    mem_width      mem_height      mem.power.total; % L3 Cache
           core1_left_x     core_bot_y     core_width     core_height     core.power.total;
           core2_left_x     core_bot_y     core_width     core_height     core.power.total;
           core3_left_x     core_bot_y     core_width     core_height     core.power.total;
           core4_left_x     core_bot_y     core_width     core_height     core.power.total];
    %blk_num is for splitting the power maps of each die
    blk_num = [8];
    
    
%%

chip_power = 4*core.power.total + mem.power.total + gpu.power.total;

package_width = 37.5e-3;
package_height = 37.5e-3;
power_therm_vec = chip_power_total;
[max_temp temp_vec] = get_stack_temperature(core.chip.num_layers,core.chip.thickness,core.wire,core.tsv,total_chip_width,total_chip_height,package_width,package_height,core.heat,simulation,map,blk_num,power_therm_vec)

    
%% Report

fprintf('Total system power consumption:\n %.4g W\n\n',chip_power)
disp('Max temperature in each tier: ')
disp(temp_vec)

    
%% Wire pitch

wire_pitch_sb_nm = [ 112.5	112.5	112.5	168.8	225	 337.6	450.1	566.5	19400 ];

figure(4)
clf
plot(wire_pitch_sb_nm,'k')
hold on
plot(core.wire.pn*1e9,'r')
xlabel('wiring layer')
ylabel('wire pitch (nm)')
grid on
ylim([0 1.2*max(core.wire.pn*1e9)])
legend('Actual','Simulated','location','nw')
fixfigs(4,3,14,12)