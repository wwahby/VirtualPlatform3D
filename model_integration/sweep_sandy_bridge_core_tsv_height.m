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
Ng_core = 86e6/4;
Ach_mm2_core = 18.5;
gate_pitch_core = 465e-9*2;
min_pitch_core = 112.5e-9;
fmax_core = 3.6e9;
w_trans = 32e-9;
Vdd_core = 1.25;

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
rent_exp_logic = 0.6;
rent_exp_mem = 0.4;
rent_exp_gpu = 0.55;

%% 
num_layers = [1 2 4 8];
layer_length = length(num_layers);

thicknesses = 1e-6:1e-6:300e-6;
num_thicknesses = length(thicknesses);

power = zeros(layer_length,num_thicknesses);
wire_power = zeros(layer_length,num_thicknesses);
temp = zeros(layer_length,num_thicknesses);
thickness = zeros(layer_length,num_thicknesses);



for nind = 1:length(num_layers)
    for thind = 1:num_thicknesses
        die_thickness = thicknesses(thind);
        num_layers_per_block = num_layers(nind);

        fprintf('====================\n')
        fprintf('==== Tiers: %d =====\n',num_layers(nind));
        fprintf('=Thickness: %d =====\n',thicknesses(thind)*1e6);
        fprintf('====================\n')

        %% define parameters

        [core.chip core.transistor core.gate core.tsv core.wire core.psn] = generate_basic_processor_settings(rent_exp_logic,num_layers_per_block,Ng_core,Ach_mm2_core,gate_pitch_core,min_pitch_core,Vdd_core,fmax_core,w_trans);

        %% Tweak wiring parameters
        core.wire.repeater_fraction = [0.3]; % 1 is default from gen_basic_proc_settings
        core.wire.routing_efficiency = [0.6]; % 0.4 is default from gen_basic_proc_settings
        core.wire.use_graphene = 0;
        simulation.force_thickness = 1;
        core.chip.thickness_nominal = die_thickness;

        %% calculate block parameters
        [core.chip core.power core.tsv core.wire core.repeater core.psn] = codesign_block(core.chip,core.tsv,core.gate,core.transistor,core.wire,heat,core.psn,simulation);

        power(nind,thind) = core.power.total;
        wire_power(nind,thind) = core.power.wiring;
        temp(nind,thind) = core.chip.temperature;
        thickness(nind,thind) = core.chip.thickness;
    end

end

    

%% Plots

figure(1)
clf
set(gcf,'DefaultAxesColorOrder',[0 0 0; 0 0 1; 1 0 0 ; 0 1 0])
hold all

for nind = 1:layer_length
    plot(thicknesses*1e6,power(nind,:))
end
xlabel('die thickness (microns)')
ylabel('Total power (W)')
fixfigs(1,3,14,12)

figure(2)
clf
set(gcf,'DefaultAxesColorOrder',[0 0 0; 0 0 1; 1 0 0 ; 0 1 0])
hold all

for nind = 1:layer_length
    plot(thicknesses*1e6,wire_power(nind,:))
end
xlabel('die thickness (microns)')
ylabel('Wiring power (W')
fixfigs(2,3,14,12)
