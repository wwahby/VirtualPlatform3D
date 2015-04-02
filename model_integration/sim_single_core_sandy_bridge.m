%% Simulation parameters
simulation.use_joyner = 0;
simulation.redo_wiring_after_repeaters = 0;
simulation.topdown_WLARI = 1; % Use topdown simultaneous WLA and RI (0 = use standard bottom-up optimal WLA, followed by one pass of RI)
simulation.skip_psn_loops = 1; % Skip PSN TSV homing for faster debug
simulation.draw_thermal_map = 0; % Plot thermal profile of each chip
simulation.print_thermal_data = 0; % Output max temp in each layer to console
simulation.separate_wiring_tiers = 1; % 1 = Each logic plane will have its own wiring tiers between it and the next logic plane
                                      % 0 = All metal layers for entire device will be routed on top of entire 3D stack

%% Logic core parameters
compression_factor = 1; % linear scaling factor. 1 = actual 32nm design, 4.57 = equivalent 7nm SB
Ng_core = 86e6/4; %86M transistors, assume 2in NAND -> /4 to get total NAND gates
Ach_mm2_core = 18.5/(compression_factor^2);
gate_pitch_core = 465e-9*2/compression_factor;
min_pitch_core = 112.5e-9/compression_factor;
fmax_core = 3.6e9;
w_trans = 32e-9/compression_factor;
Vdd_core = 1.25;

%% Thermal parameters
% %the heat transfer coefficient
% % r = 1/(hA); A is the size of top surface area
% % the cooling capability of the top heatsink; 20000, 1cm*1cm, means:
% % 0.5 W/K
% % h = q/dT - q = heat flux (W/m^2)
% heat.up = 20000;
% 
% % Bottom surface heat transfer coefficient
% % This parameter controls the area directly BELOW the bottom chip
% % If the interposer is larger than the bottom chip, heat.d controls the
% % rest of the area
% % Microfluidic heat sinks are assumed to be as large as the chip in the interposer
% heat.down = 5;  
% 
% % Heat transfer coefficient for the interposer area NOT directly underneath
% % the chip(s)
% heat.d = 5;
% 
% % Side surface heat coefficient, usually near adiabatic
% heat.side = 5;
% 
% heat.Ta = 298; % ambient temperature


r_air = 1/1.825; %K/W for a 1cm^2 HS
r_water = 1/4.63; %K/W for a 1cm^2 HS
A_hs = (1e-2)^2; % 1 cm^2

h_air = 1/(r_air*A_hs);
h_water = 1/(r_water*A_hs);
h_package = 5; % it sucks

heat.up = h_air;        % above chip
heat.down = h_package;     % directly beneath chip
heat.d = h_package;        % package, not under chip
heat.side = h_package;          % side
heat.Ta = 298; % ambient temperature

heat.interposer_thickness = 200e-6; % (m) Thickness of the interposer below the 3D stack
heat.bump_thickness = 40e-6;        % (m) Microbump thickness (between interposer and bottom chip of 3D stack)
heat.underfill_thickness = 1e-6;    % (m) Thickness of underfill material between each die in the 3D stack
heat.tim_thickness = 5e-6;          % (m) Thickness of thermal interface material between top chip in stack and heat sink

heat.material_IDs = [ 2 9 3];



%%
rent_exp_logic = 0.6;
rent_exp_mem = 0.4;
rent_exp_gpu = 0.55;

%% 
tiers = 1;
die_thickness = 10e-6;
force_thickness = 1;
rel_permittivity = 3;
frequency = fmax_core;
heat_flux_top = [ h_air ];
decap_ratio = 0.1;

t_sweep_start = cputime;

num_layers_per_block = tiers;
epsrd = rel_permittivity;

%% define parameters

[core.chip core.transistor core.gate core.tsv core.wire core.psn] = generate_basic_processor_settings(rent_exp_logic,num_layers_per_block,Ng_core,Ach_mm2_core,gate_pitch_core,min_pitch_core,Vdd_core,fmax_core,w_trans);
%core.psn.mismatch_tolerance = 0.01;
%% Tweak wiring parameters
core.wire.repeater_fraction = [0.5]; % 1 is default from gen_basic_proc_settings
core.wire.routing_efficiency = [0.5]; % 0.4 is default from gen_basic_proc_settings
core.wire.repeater_max_area_fraction = 0.3; % (-) Fraction of chip area that can be consumed by repeater/buffer gates
core.wire.repeater_via_max_area_fraction = 0.05; % (-) Fraction of routable wire area that can be consumed by vias for repeater connections

core.wire.use_graphene = 0;
simulation.force_thickness = force_thickness;
core.chip.thickness_nominal = die_thickness;
core.wire.dielectric_epsr = epsrd;
core.psn.decap_area_fraction = decap_ratio;
core.transistor.leakage_current_per_micron = 100e-9; %(A/um) % 32nm IOFF

if (die_thickness < 30e-6) % for monolithic-scale chips use thin SiO2 layer rather than underfill
    heat.interposer_thickness = 200e-6; % (m) Thickness of the interposer below the 3D stack
    heat.bump_thickness = 40e-6;        % (m) Microbump thickness (between interposer and bottom chip of 3D stack)
    heat.underfill_thickness = 0.2e-6;    % (m) Thickness of underfill material between each die in the 3D stack
    heat.tim_thickness = 5e-6;          % (m) Thickness of thermal interface material between top chip in stack and heat sink
    heat.material_IDs = [ 2 9 5];
else % for standard die stacking go ahead and use regular underfill
    heat.interposer_thickness = 200e-6; % (m) Thickness of the interposer below the 3D stack
    heat.bump_thickness = 40e-6;        % (m) Microbump thickness (between interposer and bottom chip of 3D stack)
    heat.underfill_thickness = 5e-6;    % (m) Thickness of underfill material between each die in the 3D stack
    heat.tim_thickness = 5e-6;          % (m) Thickness of thermal interface material between top chip in stack and heat sink
    heat.material_IDs = [ 2 9 3];
end

heat.up = heat_flux_top;        % above chip

%% calculate block parameters
[core.chip core.power core.tsv core.wire core.repeater core.psn] = codesign_block(core.chip,core.tsv,core.gate,core.transistor,core.wire,heat,core.psn,simulation);


t_sweep_stop = cputime;
fprintf('\nTotal time elapsed for parameter sweep: %.3g seconds\n\n',(t_sweep_stop-t_sweep_start));


%% Plots

intel_pitch = [112.5 112.5 112.5 168.8 225 337.6 450.1 566.5 19400  ]; % Actual intel data
intel_pitch_no_pow = [112.5 112.5 112.5 168.8 225 337.6 450.1 566.5];
intel_power_total = 73;

figure(1)
clf
hold on
plot(intel_pitch_no_pow,'k','linewidth',2)
plot(core.wire.pn*1e9,'b-','linewidth',2)

xlabel('Wiring tier')
ylabel('Wire Pitch (nm)')







