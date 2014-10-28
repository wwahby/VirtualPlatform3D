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
num_layers = [1 2 4 8];
thicknesses = [1e-6 100e-6];
rel_permittivities = [1 3];


layer_length = length(num_layers);
num_perms = length(rel_permittivities);
num_thicks = length(thicknesses);

power = zeros(num_thicks,layer_length,num_perms);
power_density = zeros(num_thicks,layer_length,num_perms);

wire_power = zeros(num_thicks,layer_length,num_perms);
rep_power = zeros(num_thicks,layer_length,num_perms);
temp = zeros(num_thicks,layer_length,num_perms);
thickness = zeros(num_thicks,layer_length,num_perms);
npads = zeros(num_thicks,layer_length,num_perms);
cap_density = zeros(num_thicks,layer_length,num_perms);
Ltsv_m2 = zeros(num_thicks,layer_length,num_perms);


ild_cell = cell(num_thicks,layer_length,num_perms);
psn_cell = cell(num_thicks,layer_length,num_perms);
power_cell = cell(num_thicks,layer_length,num_perms);
wire_cell = cell(num_thicks,layer_length,num_perms);
chip_cell = cell(num_thicks,layer_length,num_perms);
tsv_cell = cell(num_thicks,layer_length,num_perms);


for thind = 1:num_thicks
    for nind = 1:length(num_layers)
        for pind = 1:num_perms
            die_thickness = thicknesses(thind);
            num_layers_per_block = num_layers(nind);
            epsrd = rel_permittivities(pind);

            fprintf('====================\n')
            fprintf('==== Tiers: %d =====\n',num_layers(nind));
            fprintf('==== epsrd: %d =====\n',epsrd);
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
            core.wire.dielectric_epsr = epsrd;

            
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
            
            %% calculate block parameters
            [core.chip core.power core.tsv core.wire core.repeater core.psn] = codesign_block(core.chip,core.tsv,core.gate,core.transistor,core.wire,heat,core.psn,simulation);

            power(thind,nind,pind) = core.power.total;
            power_density(thind,nind,pind) = core.power.density;
            
            wire_power(thind,nind,pind) = core.power.wiring;
            rep_power(thind,nind,pind) = core.power.repeater;
            temp(thind,nind,pind) = core.chip.temperature;
            thickness(thind,nind,pind) = core.chip.thickness;
            npads(thind,nind,pind) = core.psn.Npads;


            ild_cell{thind,nind,pind} = core.chip.iidf;
            psn_cell{thind,nind,pind} = core.psn;
            power_cell{thind,nind,pind} = core.power;
            wire_cell{thind,nind,pind} = core.wire;
            chip_cell{thind,nind,pind} = core.chip;
            tsv_cell{thind,nind,pind} = core.tsv;

            Ltsv_m2(thind,nind,pind) = psn_cell{thind,nind,pind}.Ltsv/psn_cell{thind,nind,pind}.l_unit_cell^2;
            cap_density(thind,nind,pind) = psn_cell{thind,nind,pind}.cap_density;

        end

    end
end

%%
% for thind = 1:num_thicks
%     for pind = 1:num_perms
%         for nind = 1:length(num_layers)
%             power_density(thind,nind,pind) = power_cell{thind,nind,pind}.density;
%         end
%     end
% end


%% Plots


% We want npads(nind) for different thind and pind
% ===================
% == NPADS vs tier thickness for different ILD permittivities
% =========================================================
figure(1)
clf
%set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ])
set(gcf,'DefaultAxesColorOrder',[0 0 1; 1 0 0 ])
set(gcf,'DefaultAxesLineStyleOrder','-|--|-.')
hold all

for thind=num_thicks:-1:1
    for pind = [2 1]
    

        npads_plot = zeros(1,length(num_layers));

        for nind = 1:length(num_layers)
            npads_plot(nind) = npads(thind,nind,pind);
        end
        
        plot(num_layers,npads_plot)
    end
end
set(gca,'yscale','log')
ylim([1e1 2e3])
xlabel('Number of tiers')
ylabel('Number of Power Pads')
fixfigs(1,3,14,12)


% ===================
% == Max temp vs tier number for different ILD permittivities and
% thicknesses
% =========================================================
figure(2)
clf
%set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ])
set(gcf,'DefaultAxesColorOrder',[0 0 1; 1 0 0 ])
set(gcf,'DefaultAxesLineStyleOrder','-|--|-.')
hold all

for thind=num_thicks:-1:1
    for pind = [2 1]
    

        temp_plot = zeros(1,length(num_layers));

        for nind = 1:length(num_layers)
            temp_plot(nind) = temp(thind,nind,pind);
        end
        
        plot(num_layers,temp_plot)
        fprintf('thickness: %d \t epsrd %d\n',thicknesses(thind)*1e6,rel_permittivities(pind))
    end
end
%set(gca,'yscale','log')
xlabel('Number of tiers')
ylabel('Maximum Temperature (C)')
fixfigs(2,3,14,12)
        


% ===================
% == Power consumtion vs tier number for different ILD permittivities and
% thicknesses
% =========================================================
figure(3)
clf
%set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ])
set(gcf,'DefaultAxesColorOrder',[0 0 1; 1 0 0 ])
set(gcf,'DefaultAxesLineStyleOrder','-|--|-.')
hold all

for thind=num_thicks:-1:1
    for pind = [2 1]
        pow_plot = zeros(1,length(num_layers));

        for nind = 1:length(num_layers)
            pow_plot(nind) = power(thind,nind,pind);
        end
        
        plot(num_layers,pow_plot)
    end
end
%set(gca,'yscale','log')
xlabel('Number of tiers')
ylabel('Power (W)')
fixfigs(3,3,14,12)


% ===================
% == Power density vs tier number for different ILD permittivities and
% thicknesses
% =========================================================
figure(4)
clf
%set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1 ])
set(gcf,'DefaultAxesColorOrder',[0 0 1; 1 0 0 ])
set(gcf,'DefaultAxesLineStyleOrder','-|--|-.')
hold all

for thind=num_thicks:-1:1
    for pind = [2 1]
        pdens_plot = zeros(1,length(num_layers));

        for nind = 1:length(num_layers)
            pdens_plot(nind) = power_density(thind,nind,pind)/1e4;
        end
        
        plot(num_layers,pdens_plot)
    end
end
%set(gca,'yscale','log')
xlabel('Number of tiers')
ylabel('Power Density (W/cm^2)')
fixfigs(4,3,14,12)





