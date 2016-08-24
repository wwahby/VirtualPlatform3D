%% Simulation parameters
simulation.use_joyner = 0;
simulation.redo_wiring_after_repeaters = 0;
simulation.topdown_WLARI = 1; % Use topdown simultaneous WLA and RI (0 = use standard bottom-up optimal WLA, followed by one pass of RI)
simulation.skip_psn_loops = 0; % Skip PSN TSV homing for faster debug
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
num_layers = [1];

bulk_cu_rho = 17.2e-9;
resistivities = linspace(1/2*bulk_cu_rho,2*bulk_cu_rho,51);

% num_layers = 8;
% thicknesses = 300e-6;

layer_length = length(num_layers);
num_resistivities = length(resistivities);


power = zeros(layer_length,num_resistivities);
wire_power = zeros(layer_length,num_resistivities);
rep_power = zeros(layer_length,num_resistivities);
temp = zeros(layer_length,num_resistivities);
thickness = zeros(layer_length,num_resistivities);
npads = zeros(layer_length,num_resistivities);
num_wire_tiers = zeros(layer_length,num_resistivities);



for nind = 1:length(num_layers)
    for rind = 1:num_resistivities
        rho = resistivities(rind);
        num_layers_per_block = num_layers(nind);

        fprintf('====================\n')
        fprintf('==== Tiers: %d =====\n',num_layers(nind));
        fprintf('= Metal rho: %d =====\n',resistivities(rind)*1e9);
        fprintf('====================\n')

        %% define parameters

        [core.chip core.transistor core.gate core.tsv core.wire core.psn] = generate_basic_processor_settings(rent_exp_logic,num_layers_per_block,Ng_core,Ach_mm2_core,gate_pitch_core,min_pitch_core,Vdd_core,fmax_core,w_trans);
        core.psn.mismatch_tolerance = 0.10;
        %% Tweak wiring parameters
        core.wire.repeater_fraction = [0.3]; % 1 is default from gen_basic_proc_settings
        core.wire.routing_efficiency = [0.6]; % 0.4 is default from gen_basic_proc_settings
        core.wire.use_graphene = 0;
        
        
        core.wire.resistivity = rho;
        %simulation.force_thickness = 1;
        %core.chip.thickness_nominal = die_thickness;

        %% calculate block parameters
        [core.chip core.power core.tsv core.wire core.repeater core.psn] = codesign_block(core.chip,core.tsv,core.gate,core.transistor,core.wire,heat,core.psn,simulation);

        power(nind,rind) = core.power.total;
        wire_power(nind,rind) = core.power.wiring;
        num_wire_tiers(nind,rind) = length(core.wire.pn);
        rep_power(nind,rind) = core.power.repeater;
        temp(nind,rind) = core.chip.temperature;
        thickness(nind,rind) = core.chip.thickness;
        npads(nind,rind) = core.psn.Npads;
    end

end

    

%% Plots

figure(1)
clf
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
hold all
for nind = 1:layer_length
    plot(resistivities*1e9,power(nind,:))
end
xlabel('Metal Resistivity (Ohm m)')
ylabel('Total power (W)')
fixfigs(1,3,14,12)


figure(2)
clf
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
hold all
for nind = 1:layer_length
    plot(resistivities*1e9,num_wire_tiers(nind,:))
end
xlabel('Metal Resistivity (Ohm m)')
ylabel('Number of metal levels')
fixfigs(2,3,14,12)

figure(3)
clf
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
hold all
for nind = 1:layer_length
    plot(resistivities*1e9,wire_power(nind,:)+rep_power(nind,:))
end
xlabel('Metal Resistivity (Ohm m)')
ylabel('Total power (W)')
fixfigs(3,3,14,12)

figure(4)
clf
set(gcf,'DefaultAxesColorOrder',[0 0 1; 1 0 0; 0 0 1 ; 0 1 0])
hold all
for nind = 1:layer_length
    plot(resistivities*1e9,wire_power(nind,:))
    plot(resistivities*1e9,rep_power(nind,:))
end
xlabel('Metal Resistivity (Ohm m)')
ylabel('Power (W)')
fixfigs(4,3,14,12)

% 
% figure(2)
% clf
% set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
% hold all
% for nind = 1:layer_length
%     plot(permittivities*1e6,wire_power(nind,:))
% end
% xlabel('Layer thickness (microns)')
% ylabel('Wiring power (W)')
% set(gca,'xscale','log')
% xlim([1 300])
% fixfigs(2,3,14,12)
% 
% figure(3)
% clf
% %set(gcf,'DefaultAxesColorOrder',[1 0 0; 0 0 1 ; 0 1 0])
% set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
% hold all
% for nind = 1:layer_length
%     plot(permittivities*1e6,npads(nind,:)*2)
% end
% xlabel('Layer thickness (microns)')
% ylabel('Number of power delivery TSVs')
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% xlim([1 300])
% fixfigs(3,3,14,12)