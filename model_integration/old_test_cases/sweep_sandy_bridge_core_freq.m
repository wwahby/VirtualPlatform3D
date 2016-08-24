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

heat.up = h_air;           % Above the chip
heat.down = h_package;     % Portion of the package that is directly under the chip
heat.d = h_package;        % Portion of the package that is not directly under chip
heat.side = h_package;     % Side of the chip
heat.Ta = 298;             % (K) ambient temperature

%%
rent_exp_logic = 0.6;
rent_exp_mem = 0.4;
rent_exp_gpu = 0.55;

%% 
num_layers = [1 2 4 8];

frequencies = logspace(8,10,21);
frequencies = [frequencies 20e9 30e9];

%thicknesses = [1e-6];% 10e-6 50e-6 100e-6 200e-6 300e-6];

% num_layers = 8;
% thicknesses = 300e-6;

layer_length = length(num_layers);
num_freqs = length(frequencies);


power = zeros(layer_length,num_freqs);
wire_power = zeros(layer_length,num_freqs);
rep_power = zeros(layer_length,num_freqs);
temp = zeros(layer_length,num_freqs);
thickness = zeros(layer_length,num_freqs);
npads = zeros(layer_length,num_freqs);
ild_mat = cell(layer_length,num_freqs);


die_thickness = 1e-6; % (m)
for nind = 1:length(num_layers)
    for freqind = 1:num_freqs

        fmax_core = frequencies(freqind);
        num_layers_per_block = num_layers(nind);

        fprintf('========================\n')
        fprintf('==== Tiers: %d =========\n',num_layers(nind));
        fprintf('=Frequency: %d GHz =====\n',frequencies(freqind)/1e9);
        fprintf('========================\n')

        %% define parameters

        [core.chip core.transistor core.gate core.tsv core.wire core.psn] = generate_basic_processor_settings(rent_exp_logic,num_layers_per_block,Ng_core,Ach_mm2_core,gate_pitch_core,min_pitch_core,Vdd_core,fmax_core,w_trans);
        core.psn.mismatch_tolerance = 0.10;
        %% Tweak wiring parameters
        core.wire.repeater_fraction = [0.3]; % 1 is default from gen_basic_proc_settings
        core.wire.routing_efficiency = [0.6]; % 0.4 is default from gen_basic_proc_settings
        core.wire.use_graphene = 0;
        simulation.force_thickness = 1;
        core.chip.thickness_nominal = die_thickness;

        %% calculate block parameters
        [core.chip core.power core.tsv core.wire core.repeater core.psn] = codesign_block(core.chip,core.tsv,core.gate,core.transistor,core.wire,heat,core.psn,simulation);

        power(nind,freqind) = core.power.total;
        wire_power(nind,freqind) = core.power.wiring;
        rep_power(nind,freqind) = core.power.repeater;
        temp(nind,freqind) = core.chip.temperature;
        thickness(nind,freqind) = core.chip.thickness;
        npads(nind,freqind) = core.psn.Npads;
        ild_mat{nind,freqind} = core.chip.iidf;

    end

end

%% gather all the results in case we want to save them for later use
results.power = power;
results.wire_power = wire_power;
results.rep_power = rep_power;
results.temp = temp;
results.thickness = thickness;
results.npads = npads;
results.ild_mat = ild_mat;
results.frequencies = frequencies;

%% Find the maximum frequency that the device can run to stay below 90C
Tmax = 90;
fmax_temp_limited = zeros(1,layer_length);

for nind =  1:layer_length
    Tind = find( (temp(nind,:) <= Tmax), 1, 'last');
    fmax_temp_limited(nind) = frequencies(Tind);
end

%% Plots

% Total power vs frequency
figure(1)
clf
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
hold all
for nind = 1:layer_length
    plot(frequencies/1e9,power(nind,:))
end
xlabel('Clock frequency (GHz)')
ylabel('Total power (W)')
set(gca,'xscale','log')
%set(gca,'yscale','log')
fixfigs(1,3,14,12)

% Wire power (ignores repeaters) vs frequency
figure(2)
clf
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
hold all
for nind = 1:layer_length
    plot(frequencies/1e9,wire_power(nind,:))
end
xlabel('Clock frequency (GHz)')
ylabel('Wiring power (W)')
set(gca,'xscale','log')
xlim([0.1 10])
fixfigs(2,3,14,12)

% Power pads vs frequency
figure(3)
clf
%set(gcf,'DefaultAxesColorOrder',[1 0 0; 0 0 1 ; 0 1 0])
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
hold all
for nind = 1:layer_length
    plot(frequencies/1e9,npads(nind,:)*2)
end
xlabel('Clock frequency (GHz)')
ylabel('Number of power delivery TSVs')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([0.1 10])
fixfigs(3,3,14,12)

% On-chip communication power vs frequency
figure(4)
clf
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
hold all
for nind = 1:layer_length
    plot(frequencies/1e9,wire_power(nind,:) + rep_power(nind,:))
end
xlabel('Clock frequency (GHz)')
ylabel('On-chip communication power (W)')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([0.1 10])
fixfigs(4,3,14,12)

% Maximum temperature vs frequency
figure(5)
clf
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
hold all
for nind = 1:layer_length
    plot(frequencies/1e9,temp(nind,:))
end
xlabel('Clock frequency (GHz)')
ylabel('Maximum temperature (C)')
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([0.1 10])
fixfigs(5,3,14,12)

% On-chip communication power as a fraction of total power vs frequency
comm_power = wire_power + rep_power;
comm_power_frac = comm_power./power;
figure(8)
clf
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
hold all
for nind = 1:layer_length
    plot(frequencies/1e9,comm_power_frac(nind,:))
end
xlabel('Clock frequency (GHz)')
ylabel('On-chip Communication Power Fraction')
set(gca,'xscale','log')
%set(gca,'yscale','log')
xlim([0.1 10])
ylim([0 1])
fixfigs(8,3,14,12)
%%

figure(6)
clf
set(gcf,'DefaultAxesColorOrder',[ 1 0 0; 1 0 1; 0 1 0 ; 0 0 1])
hold all
nind = layer_length;

for freqind = fliplr([1 3 4 5])
    plot(ild_mat{nind,freqind})
end
plot(ild_mat{1,1},'k')
set(gca,'yscale','log')
ylim([1e-1 1e3])
xlabel('Wire length (GP)')
ylabel('Expected number of on-chip wires')
fixfigs(6,5,18,16)


%% 
figure(7)
clf
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
hold all
%for nind = 1:layer_length
    plot(num_layers,power(:,1))
    plot(num_layers,wire_power(:,1)+rep_power(:,1))
    plot(num_layers, power(:,1) - (wire_power(:,1)+rep_power(:,1)) )
%end
xlabel('Number of tiers')
ylabel('Power (W)')
fixfigs(7,3,14,12)