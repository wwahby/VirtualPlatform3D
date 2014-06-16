close all

%% Simulation parameters
simulation.use_joyner = 0;
simulation.redo_wiring_after_repeaters = 0;
simulation.topdown_WLARI = 1; % Use topdown simultaneous WLA and RI (0 = use standard bottom-up optimal WLA, followed by one pass of RI)
simulation.skip_psn_loops = 1; % Skip PSN TSV homing for faster debug
simulation.draw_thermal_map = 0; % Plot thermal profile of each chip
simulation.print_thermal_data = 0; % Output max temp in each layer to console
simulation.separate_wiring_tiers = 1; % 1 = calculate wire pitch for wiring tiers between EACH logic plane

%% Logic core parameters

% Ng_core = 86e6/4;
% Ach_mm2_core = 18.5;
% gate_pitch_core = 465e-9*2;
% min_pitch_core = 112.5e-9;
% fmax_core = 3.5e9; % 3.5GHz normal, 3.9GHz turbo
% w_trans = 32e-9;
% Vdd_core = 1.25;

scale_factor = 4;

Ng_core = 86e6/4;
Ach_mm2_core = 18.5/scale_factor^2;
gate_pitch_core = 465e-9*2/scale_factor;
%min_pitch_core = 112.5e-9/scale_factor;
min_pitch_core = 5e-9;
fmax_core = 3.9e9;
w_trans = 32e-9/scale_factor;
Vdd_core = 1.25;
use_graphene = 0;

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

%% define parameters

[core.chip core.transistor core.gate core.tsv core.wire core.psn] = generate_basic_processor_settings(rent_exp_logic,num_layers_per_block,Ng_core,Ach_mm2_core,gate_pitch_core,min_pitch_core,Vdd_core,fmax_core,w_trans);

%% Tweak wiring parameters
core.wire.repeater_fraction = [0.3]; % 1 is default from gen_basic_proc_settings
core.wire.routing_efficiency = [0.6]; % 0.4 is default from gen_basic_proc_settings

%% calculate block parameters
core.wire.use_graphene = 0;
[cu_core.chip cu_core.power cu_core.tsv cu_core.wire cu_core.repeater cu_core.psn] = codesign_block(core.chip,core.tsv,core.gate,core.transistor,core.wire,heat,core.psn,simulation);

core.wire.use_graphene = 1;
[gnr_core.chip gnr_core.power gnr_core.tsv gnr_core.wire gnr_core.repeater gnr_core.psn] = codesign_block(core.chip,core.tsv,core.gate,core.transistor,core.wire,heat,core.psn,simulation);

    
%% Report

fprintf('Cu core power consumption:\n %.4g W\n\n',cu_core.power.total)
fprintf('GNR core power consumption:\n %.4g W\n\n',gnr_core.power.total)

cu_core.wire.material_vec
gnr_core.wire.material_vec

    
%% Wire pitch

wire_pitch_sb_nm = [ 112.5	112.5	112.5	168.8	225	 337.6	450.1	566.5	19400 ];

figure(4)
clf
plot(cu_core.wire.pn*1e9,'r')
hold on
plot(gnr_core.wire.pn*1e9,'b')
plot(gnr_core.wire.pn_cu*1e9,'b--')
xlabel('wiring layer')
ylabel('wire pitch (nm)')
grid on
ylim([0 150])
fixfigs(4,3,14,12)

figure(5)
clf
semilogy(cu_core.wire.Ln,'r')
hold on
semilogy(gnr_core.wire.Ln,'b')
xlabel('wiring layer')
ylabel('Longest wire routed (GP)')
fixfigs(5,3,14,12)
