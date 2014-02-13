%% Simulation parameters
simulation.use_joyner = 0;
simulation.redo_wiring_after_repeaters = 0;
simulation.topdown_WLARI = 1; % Use topdown simultaneous WLA and RI (0 = use standard bottom-up optimal WLA, followed by one pass of RI)
simulation.skip_psn_loops = 1; % Skip PSN TSV homing for faster debug
simulation.draw_thermal_map = 0; % Plot thermal profile of each chip
simulation.print_thermal_data = 0; % Output max temp in each layer to console

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
Ach_mm2_mem = 62;
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
core_height_rel = 2.5;
core_width = sqrt(Ach_mm2_core*1e-6/core_height_rel);
core_height = core_width*core_height_rel;

mem_width_rel = 4;
mem_height = sqrt(Ach_mm2_mem*1e-6/mem_width_rel);
mem_width = mem_height*mem_width_rel;

chip_width = mem_width;
chip_height = mem_height + core_height;

    %power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....
    map = [0     0     mem_width    mem_height     mem.power.total;
           0     mem_height  core_width    core_height     core.power.total;
           core_width     mem_height  core_width    core_height     core.power.total;
           2*core_width     mem_height  core_width    core_height     core.power.total;
           3*core_width     mem_height  core_width    core_height     core.power.total];
    %blk_num is for splitting the power maps of each die
    blk_num = [5];
    
    
%%

chip_power = 4*core.power.total + mem.power.total;
simulation.draw_thermal_map = 1; % Plot thermal profile of each chip
simulation.print_thermal_data = 1; % Output max temp in each layer to console

%% Thermal module -- Find actual system temperature

%%%%%%%%%%%%%%%%%%%%%%%geometry information of the chip%%%%%%%%%%%%%%%%%%%
    die.N = core.chip.num_layers;
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
    chip_therm.Xsize = chip_width; %x dimension of chip
    chip_therm.Ysize = chip_height; %y dimension of chip
    chip_therm.Xgrid = chip_therm.Xsize/grid_factor; %x grid size of chip
    chip_therm.Ygrid = chip_therm.Ysize/grid_factor; %y grid size of chip
    
    %for the interposer dimension
    pack.Xsize = 3.5*chip_width;
    pack.Ysize = 3.5*chip_height;
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
    
    power_therm_vec = ones(1,core.chip.num_layers)*chip_power/core.chip.num_layers;  %power dissipation of each die
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
