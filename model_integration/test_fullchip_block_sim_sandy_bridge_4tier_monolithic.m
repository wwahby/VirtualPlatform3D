%% Simulation parameters
simulation.use_joyner = 0;
simulation.redo_wiring_after_repeaters = 0;
simulation.topdown_WLARI = 1; % Use topdown simultaneous WLA and RI (0 = use standard bottom-up optimal WLA, followed by one pass of RI)
simulation.skip_psn_loops = 1; % Skip PSN TSV homing for faster debug
simulation.draw_thermal_map = 0; % Plot thermal profile of each chip
simulation.print_thermal_data = 0; % Output max temp in each layer to console
simulation.separate_wiring_tiers = 0; % 1 = Each logic plane will have its own wiring tiers between it and the next logic plane
                                      % 0 = All metal layers for entire device will be routed on top of entire 3D stack

%% Logic core parameters

% Ng_core = 86e6/4;
% Ach_mm2_core = 18.5;
% gate_pitch_core = 932e-9;
% min_pitch_core = 112.5e-9;
% fmax_core = 3.9e9; % 3.5GHz normal, 3.9GHz turbo
% w_trans = 32e-9;
% Vdd_core = 1.25;

Ng_core = 86e6/4;
Ach_mm2_core = 18.5;
gate_pitch_core = 465e-9*2;
min_pitch_core = 112.5e-9;
fmax_core = 3.6e9;
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
fmax_gpu = 1.35e9;
Vdd_gpu = 1.25;

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
num_layers_per_block = 4;

rent_exp_logic = 0.6;
rent_exp_mem = 0.4;
rent_exp_gpu = 0.55;

%% define parameters

[core.chip core.transistor core.gate core.tsv core.wire core.psn] = generate_basic_processor_settings(rent_exp_logic,num_layers_per_block,Ng_core,Ach_mm2_core,gate_pitch_core,min_pitch_core,Vdd_core,fmax_core,w_trans);
[mem.chip mem.transistor mem.gate mem.tsv mem.wire mem.psn] = generate_basic_processor_settings(rent_exp_mem,num_layers_per_block,Ng_mem,Ach_mm2_mem,gate_pitch_mem,min_pitch_mem,Vdd_mem,fmax_mem,w_trans);
[gpu.chip gpu.transistor gpu.gate gpu.tsv gpu.wire gpu.psn] = generate_basic_processor_settings(rent_exp_gpu,num_layers_per_block,Ng_gpu,Ach_mm2_gpu,gate_pitch_gpu,min_pitch_gpu,Vdd_gpu,fmax_gpu,w_trans);

%% Tweak wiring parameters
core.wire.repeater_fraction = [0.3]; % 1 is default from gen_basic_proc_settings
core.wire.routing_efficiency = [0.6]; % 0.4 is default from gen_basic_proc_settings
core.wire.use_graphene = 0;

gpu.wire.repeater_fraction = core.wire.repeater_fraction;
gpu.wire.routing_efficiency = core.wire.routing_efficiency;
gpu.wire.use_graphene = core.wire.use_graphene;

mem.wire.repeater_fraction = core.wire.repeater_fraction;
mem.wire.routing_efficiency = core.wire.routing_efficiency;
mem.wire.use_graphene = core.wire.use_graphene;

%% calculate block parameters
[core.chip core.power core.tsv core.wire core.repeater core.psn] = codesign_block(core.chip,core.tsv,core.gate,core.transistor,core.wire,heat,core.psn,simulation);
[mem.chip mem.power mem.tsv mem.wire mem.repeater mem.psn] = codesign_block(mem.chip,mem.tsv,mem.gate,mem.transistor,mem.wire,heat,mem.psn,simulation);
[gpu.chip gpu.power gpu.tsv gpu.wire gpu.repeater gpu.psn] = codesign_block(gpu.chip,gpu.tsv,gpu.gate,gpu.transistor,gpu.wire,heat,gpu.psn,simulation);


%% Define block dimensions (m)
chip_width = 21.08e-3;
chip_height = 10.25e-3;

core_width = 3.302e-3;
core_height = 5.65e-3;
core_area = core_width*core_height;
core_width = sqrt(core_area);
core_height = core_width;

gpu_width = 4.89e-3;
gpu_height = 9.12e-3;
gpu_area = gpu_width * gpu_height;
gpu_height = 2*core_height;
gpu_width = gpu_area/gpu_height;

mem_width = 13.24e-3;
mem_height = 3.48e-3;
mem_area = mem_width * mem_height;
mem_width = 2*core_width + gpu_width;
mem_height = mem_area/mem_width;

mem_io_width = 13.49e-3;
mem_io_height = 1.16e-3;
mem_io_area = mem_io_width * mem_io_height;
mem_io_width = mem_width;
mem_io_height = mem_io_area/mem_io_width;

sys_agent_width = 3.01e-3;
sys_agent_height = 9.11e-3;
sys_agent_area = sys_agent_width * sys_agent_height;
sys_agent_width = mem_width;
sys_agent_height = sys_agent_area/sys_agent_width;

core_bot_y_12 = 0;
core_bot_y_34 = core_height;
core1_left_x = gpu_width;
core2_left_x = core1_left_x + core_width;
core3_left_x = gpu_width;
core4_left_x = core3_left_x + core_width;

gpu_left_x = 0;
gpu_bot_y = 0;

mem_left_x = 0;
mem_bot_y = mem_io_height;

mem_io_left_x = 0;
mem_io_bot_y = 0;

sys_agent_left_x = 0;
sys_agent_bot_y = mem_io_height + mem_height;


%% Set up block map
    %power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....
    map =  [gpu_left_x      gpu_bot_y       gpu_width      gpu_height      gpu.power.total; % GPU
           core1_left_x     core_bot_y_12     core_width     core_height     core.power.total;
           core2_left_x     core_bot_y_12     core_width     core_height     core.power.total;
           core3_left_x     core_bot_y_34     core_width     core_height     core.power.total;
           core4_left_x     core_bot_y_34     core_width     core_height     core.power.total;
           mem_io_left_x   mem_io_bot_y    mem_io_width   mem_io_height     0; % memory IO
           sys_agent_left_x   sys_agent_bot_y    sys_agent_width   sys_agent_height     0; % system agent + memory controller
           mem_left_x     mem_io_height    mem_width      mem_height      mem.power.total ];
    %blk_num is for splitting the power maps of each die
    blk_num = [5 3];
    
    
%%

chip_power = 4*core.power.total + mem.power.total + gpu.power.total;
power_therm_vec = [(4*core.power.total + gpu.power.total) mem.power.total];  %power dissipation of each die


simulation.draw_thermal_map = 1; % Plot thermal profile of each chip
simulation.print_thermal_data = 1; % Output max temp in each layer to console

%% Thermal module -- Find actual system temperature

%%%%%%%%%%%%%%%%%%%%%%%geometry information of the chip%%%%%%%%%%%%%%%%%%%
    die.N = 2;
    die.model = 3; % for each die, how many layers we model
    
    % flip chip package; order:
    %heatsink->TIM->CHIP_BULK1->METAL->BONDING->CHIP_BULK2-> ...
    %               CHIP_BULKN->METAL->MICRO-BUMPS->INTERPOSER    
    thick.bump = 40e-6;  %micro-bump thickness; between second die and interposer
    thick.tim = 5e-6; %tim thickness; between the chip and heatsink
    thick.under = 5e-6; %underfill bonding thickness; between two dies
    thick.inter = 200e-6; %interposer thickness
    thick.die = core.chip.thickness; %die thickness
    thick.ild = sum(core.wire.pn); %metal layer thickness

    grid_factor = 50;
    chip_therm.Xsize = 2*core_width + gpu_width; %x dimension of chip
    chip_therm.Ysize = 2*core_height; %y dimension of chip
    chip_therm.Xgrid = chip_therm.Xsize/grid_factor; %x grid size of chip
    chip_therm.Ygrid = chip_therm.Ysize/grid_factor; %y grid size of chip
    
    %for the interposer dimension
    pack.Xsize = 37.5e-3;
    pack.Ysize = 37.5e-3;
    pack.Xgrid = pack.Xsize/grid_factor;
    pack.Ygrid = pack.Ysize/grid_factor;
    
    %assumed TSV starting from top metal layer of a top die
    %            to the first metal layer of a bottom die
    %Thus the TSV passes through bonding layer & bulk of a botom die
    %TSV geometry 
    tsv_therm.d = core.tsv.width_m; % tsv diameter including the liner thickness
    tsv_therm.liner = core.tsv.width_m/10; %liner thickness
    tsv_therm.px = core.tsv.pitch_m; %x direction pitch
    tsv_therm.py = core.tsv.pitch_m; %y direction pitch
    tsv_therm.Nx = round(sqrt(core.tsv.num)); %x direction number
    tsv_therm.Ny = round(sqrt(core.tsv.num)); %y direction number

    %Bump geometry
    bump.d = 20e-6;
    bump.px = 100e-6;
    bump.py = 100e-6;
    bump.Nx = 40;
    bump.Ny = 40;
    
    portion = 0.5; %the metal portion in the ILD layers
    %This is used for equivalent thermal resistance calculation of metal
    %layer
   
    % from top to bottom; unit: watt
    
    granularity = 20;
    % thermal map, number of color used
    
    draw = simulation.draw_thermal_map;
    drawP = simulation.draw_thermal_map;
    % whether to draw the thermal map; 1 yes; 0 no
    
    displayT = simulation.print_thermal_data;
    % print the temperature information


    chip.temperature_vec = thermal.ThermSim( die, thick, chip_therm, pack, ...
              tsv_therm, bump, portion, power_therm_vec, ...
              map, blk_num, granularity, draw, drawP, heat, displayT);
          
    chip.temperature = max(chip.temperature_vec);
    
%% Report

fprintf('Total system power consumption:\n %.4g W\n\n',chip_power)
disp('Max temperature in each tier: ')
disp(chip.temperature_vec)

    
%% Wire pitch

wire_pitch_sb_nm = [ 112.5	112.5	112.5	168.8	225	 337.6	450.1	566.5	19400 ];

figure(9)
clf
plot(wire_pitch_sb_nm,'k')
hold on
plot(core.wire.pn*1e9,'r')
xlabel('wiring layer')
ylabel('wire pitch (nm)')
grid on
ylim([0 700])
fixfigs(9,3,14,12)