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

heat.up = h_air;        % above chip
heat.down = h_package;     % directly beneath chip
heat.d = h_package;        % package, not under chip
heat.side = h_package;          % side
heat.Ta = 298; % ambient temperature

%%
rent_exp_logic = 0.6;
rent_exp_mem = 0.4;
rent_exp_gpu = 0.55;

%% 
num_layers = [1 2 4 8];
decap_ratios = logspace(-3,-1,21);


layer_length = length(num_layers);
num_decaps = length(decap_ratios);

power = zeros(layer_length,num_decaps);
wire_power = zeros(layer_length,num_decaps);
rep_power = zeros(layer_length,num_decaps);
temp = zeros(layer_length,num_decaps);
thickness = zeros(layer_length,num_decaps);
npads = zeros(layer_length,num_decaps);
cap_density = zeros(layer_length,num_decaps);
Ltsv_m2 = zeros(layer_length,num_decaps);


ild_cell = cell(layer_length,num_decaps);
psn_cell = cell(layer_length,num_decaps);
power_cell = cell(layer_length,num_decaps);
wire_cell = cell(layer_length,num_decaps);
chip_cell = cell(layer_length,num_decaps);
tsv_cell = cell(layer_length,num_decaps);


for nind = 1:length(num_layers)
    for dind = 1:num_decaps
        die_thickness = 100e-6;
        num_layers_per_block = num_layers(nind);
        decap_ratio = decap_ratios(dind);

        fprintf('====================\n')
        fprintf('==== Tiers: %d =====\n',num_layers(nind));
        fprintf('==== Decap: %d =====\n',decap_ratio);
        fprintf('====================\n')

        %% define parameters

        [core.chip core.transistor core.gate core.tsv core.wire core.psn] = generate_basic_processor_settings(rent_exp_logic,num_layers_per_block,Ng_core,Ach_mm2_core,gate_pitch_core,min_pitch_core,Vdd_core,fmax_core,w_trans);
        %core.psn.mismatch_tolerance = 0.01;
        %% Tweak wiring parameters
        core.wire.repeater_fraction = [0.3]; % 1 is default from gen_basic_proc_settings
        core.wire.routing_efficiency = [0.6]; % 0.4 is default from gen_basic_proc_settings
        core.wire.use_graphene = 0;
        simulation.force_thickness = 1;
        core.chip.thickness_nominal = die_thickness;
        core.psn.decap_area_fraction = decap_ratio;

        %% calculate block parameters
        [core.chip core.power core.tsv core.wire core.repeater core.psn] = codesign_block(core.chip,core.tsv,core.gate,core.transistor,core.wire,heat,core.psn,simulation);

        power(nind,dind) = core.power.total;
        wire_power(nind,dind) = core.power.wiring;
        rep_power(nind,dind) = core.power.repeater;
        temp(nind,dind) = core.chip.temperature;
        thickness(nind,dind) = core.chip.thickness;
        npads(nind,dind) = core.psn.Npads;
        
        
        ild_cell{nind,dind} = core.chip.iidf;
        psn_cell{nind,dind} = core.psn;
        power_cell{nind,dind} = core.power;
        wire_cell{nind,dind} = core.wire;
        chip_cell{nind,dind} = core.chip;
        tsv_cell{nind,dind} = core.tsv;
        
        Ltsv_m2(nind,dind) = psn_cell{nind,dind}.Ltsv/psn_cell{nind,dind}.l_unit_cell^2;
        cap_density(nind,dind) = psn_cell{nind,dind}.cap_density;
        
    end

end


%%

% cind = 2;
% Jch = power_cell{cind}.density/chip_cell{cind}.Vdd;

fprintf('\n\n')
fprintf('Tiers \t Decap \t Npads \t Cd_m2 \t L_m2\n')
for nind = 1:length(num_layers)
    for dind = 1:num_decaps
        fprintf('%d \t\t %.3g \t %.3d \t %.3d \t %.3d\n',num_layers(nind),decap_ratios(dind),psn_cell{nind,dind}.Npads,psn_cell{nind,dind}.cap_density,psn_cell{nind,dind}.Ltsv/psn_cell{nind,dind}.l_unit_cell^2)
    end
end
        
    

%% Decap plot: Npads vs C
% ===========================
figure(1)
clf
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
hold all
for nind = 1:length(num_layers)
    npads_n = zeros(1,num_decaps);
    for dind = 1:num_decaps
        npads_n(dind) = psn_cell{nind,dind}.Npads;
    end
    plot(cap_density(nind,:),npads_n)
end
set(gca,'yscale','log')
set(gca,'xscale','log')
xlim([min(cap_density(1,:)) max(cap_density(1,:))])
xlabel('Capacitance Density (F/m^2)')
ylabel('Number of Power Pads Required')
fixfigs(1,3,14,12)

% ===========================
% Decap: Npads vs Decap ratio
% ===========================
figure(2)
clf
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
hold all
for nind = 1:length(num_layers)
    npads_n = zeros(1,dind);
    for dind = 1:num_decaps
        npads_n(dind) = psn_cell{nind,dind}.Npads;
    end
    plot(decap_ratios,npads_n)
end
set(gca,'yscale','log')
set(gca,'xscale','log')
xlim([1e-3 max(decap_ratios)])
%ylim([1e1 1e4])
xlabel('Decap ratio')
ylabel('Number of Power Pads Required')
fixfigs(2,3,14,12)



% ===========================
% Decap Plot: Cdens vs Ldens
% ===========================
figure(3)
clf
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
set(gcf,'DefaultAxesLineStyleOrder','-')
hold all
for nind = 1:length(num_layers)
    plot(decap_ratios,cap_density(nind,:))
end

set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
set(gcf,'DefaultAxesLineStyleOrder','--')

for nind = 1:length(num_layers)
    Ldens_n = zeros(1,dind);
    for dind = 1:num_decaps
        
        
        if(nind == 1)
            psn_cell{nind,dind}.Ltsv = 0; % 2D incorrectly reports infinite TSV inductance
        end
        
        Ltsv = psn_cell{nind,dind}.Ltsv;
        Lpkg = psn_cell{nind,dind}.package_inductance;
        Ltot = Ltsv + Lpkg;
        
        Ldens_n(dind) = Ltot/psn_cell{nind,dind}.l_unit_cell^2;
    end
    plot(decap_ratios,Ldens_n,'--')
end
set(gca,'yscale','log')
xlabel('Decap ratio')
ylabel('Reactance Density')
fixfigs(3,3,14,12)

% ===========================
% Decap Plot: C vs L
% ===========================
figure(4)
clf
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
set(gcf,'DefaultAxesLineStyleOrder','-')
hold all
for nind = 1:length(num_layers)
    C_n = zeros(1,dind);
    for dind = 1:num_decaps
        C_n(dind) = cap_density(nind,dind) *  chip_cell{nind,dind}.area_total;
    end
    plot(decap_ratios,C_n)
end

set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
set(gcf,'DefaultAxesLineStyleOrder','--')

for nind = 1:length(num_layers)
    L_n = zeros(1,dind);
    for dind = 1:num_decaps
        
        if(nind == 1)
            psn_cell{nind,dind}.Ltsv = 0; % 2D incorrectly reports infinite TSV inductance
        end
        
        Ltsv = psn_cell{nind,dind}.Ltsv;
        Lpkg = psn_cell{nind,dind}.package_inductance;
        Ltot = Ltsv + Lpkg;
        
        Ldens_n(dind) = Ltot/psn_cell{nind,dind}.l_unit_cell^2;
        L_n(nind) = Ldens_n(nind) * chip_cell{nind,dind}.area_total;
    end
    plot(decap_ratios,L_n,'--')
end
set(gca,'yscale','log')
xlabel('Decap ratio')
ylabel('Absolute Reactance')
fixfigs(4,3,14,12)

% ===========================
% Decap Plot: C/ L
% ===========================
figure(5)
clf
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
set(gcf,'DefaultAxesLineStyleOrder','-')
hold all
for nind = 1:length(num_layers)
    C_n = zeros(1,dind);
    L_n = zeros(1,dind);
    for dind = 1:num_decaps
        if(nind == 1)
            psn_cell{nind,dind}.Ltsv = 0; % 2D incorrectly reports infinite TSV inductance
        end
        
        Ltsv = psn_cell{nind,dind}.Ltsv;
        Lpkg = psn_cell{nind,dind}.package_inductance;
        Ltot = Ltsv + Lpkg;
        C_n(dind) = cap_density(nind,dind) *  chip_cell{nind,dind}.area_total;
        Ldens_n(dind) = Ltot/psn_cell{nind,dind}.l_unit_cell^2;
        L_n(nind) = Ldens_n(nind) * chip_cell{nind,dind}.area_total;
        chip_cell{nind,dind}.area_total;
    end
    plot(decap_ratios,Ldens_n./cap_density(nind,:))
end

%plot(decap_ratios,ones(1,length(decap_ratios)),'Color',[0.2 0.2 0.2])

xlabel('Decap ratio')
ylabel('L/C')
set(gca,'yscale','log')
set(gca,'xscale','log')
xlim([min(decap_ratios) max(decap_ratios)])
fixfigs(5,3,14,12) 

% ===========================
% Fill factors
% ===========================
figure(6)
clf
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
set(gcf,'DefaultAxesLineStyleOrder','-')
hold all
for nind = 1:length(num_layers)
    bot_ff = zeros(1,num_decaps);
    top_ff = zeros(1,num_decaps);
    for dind = 1:num_decaps
        Awires = wire_cell{nind,dind}.wire_area;
        Avias = wire_cell{nind,dind}.via_area;
        Aused = Awires + Avias;
        Aavail = wire_cell{nind,dind}.area_per_layer;
        %Aavail_last = wire_cell{nind,dind}.routing_efficiency(end)*Aavail(1);
        ff_last = Aused(end)/Aavail(1);
        ff = Aused./Aavail;
        bot_ff(dind) = ff(1);
        top_ff(dind) = ff_last;
    end
    
    plot(decap_ratios,1-bot_ff,'-');
    %plot(decap_ratios,top_ff,'--');
end

%set(gcf,'DefaultAxesLineStyleOrder','--')
for nind = 1:length(num_layers)
    bot_ff = zeros(1,num_decaps);
    top_ff = zeros(1,num_decaps);
    for dind = 1:num_decaps
        Awires = wire_cell{nind,dind}.wire_area;
        Avias = wire_cell{nind,dind}.via_area;
        Aused = Awires + Avias;
        Aavail = wire_cell{nind,dind}.area_per_layer;
        %Aavail_last = wire_cell{nind,dind}.routing_efficiency(end)*Aavail(1);
        ff_last = Aused(end)/Aavail(1);
        ff = Aused./Aavail;
        bot_ff(dind) = ff(1);
        top_ff(dind) = ff_last;
    end
    
    %plot(decap_ratios,bot_ff);
    plot(decap_ratios,top_ff,'--');
end
xlabel('Decap ratio')
ylabel('Wiring Fill Factor')
ylim([1e-3 1e0])
set(gca,'yscale','log')
fixfigs(6,3,14,12)

% ===========================
% TSVs
% ===========================
figure(7)
clf
hold all
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
set(gcf,'DefaultAxesLineStyleOrder','-')
for nind = 1:length(num_layers)
    wtsv_vec = zeros(1,dind);
    for dind = 1:num_decaps
        wtsv_vec(dind) = tsv_cell{nind,dind}.width_gp;
    end
    plot(decap_ratios,wtsv_vec);
end
xlabel('Decap ratio')
ylabel('TSV Width (GP)')
%ylim([1e-3 1e0])
%set(gca,'yscale','log')
fixfigs(7,3,14,12)


% ===========================
% Number of wiring tiers
% ===========================
figure(8)
clf
hold all
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
set(gcf,'DefaultAxesLineStyleOrder','-')
for nind = 1:length(num_layers)
    num_tiers = zeros(1,num_decaps);
    for dind = 1:num_decaps
        num_tiers(dind) = length(wire_cell{nind,dind}.pn);
    end
    plot(decap_ratios,num_tiers);
end
xlabel('Decap ratio')
ylabel('Number of wiring tiers')
fixfigs(8,3,14,12)



% ===========================
% PSN Noise
% ===========================
figure(9)
clf
hold all
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
set(gcf,'DefaultAxesLineStyleOrder','-')
psn_target = psn_cell{1,1}.noise_target;
psn_tol = psn_cell{1,1}.mismatch_tolerance;

for nind = 1:length(num_layers)
    psn_max_vec = zeros(1,num_decaps);
    for dind = 1:num_decaps
        psn_max_vec(dind) = psn_cell{nind,dind}.noise;
    end
    
    plot(decap_ratios,psn_max_vec)
end
set(gca,'xscale','log')
xlabel('Decap ratio')
ylabel('PSN Noise (V)')
fixfigs(9,3,14,12)



% ===========================
% PSN Rel Error
% ===========================
figure(10)
clf
hold all
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ; 0 1 0])
set(gcf,'DefaultAxesLineStyleOrder','-')
psn_target = psn_cell{1,1}.noise_target;
psn_tol = psn_cell{1,1}.mismatch_tolerance;

for nind = 1:length(num_layers)
    psn_max_vec = zeros(1,num_decaps);
    for dind = 1:num_decaps
        psn_max_vec(dind) = psn_cell{nind,dind}.noise;
    end
    
    err = psn_target - psn_max_vec;
    rel_err = err/psn_target;
    norm_err = abs(rel_err);
    plot(decap_ratios,norm_err)
end
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('Decap ratio')
ylabel('Relative error')
fixfigs(10,3,14,12)
        





