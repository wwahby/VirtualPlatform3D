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
tiers = [1 2 4 8];
thicknesses = [1e-6];
force_thickness = 1;
rel_permittivities = 1:4;
frequencies = 1e9;
heat_fluxes = [ h_air ];
decap_ratios = [0.1];
wire_resistivities = [17.2e-9];


num_stacks = length(tiers);
num_perms = length(rel_permittivities);
num_thicks = length(thicknesses);
num_freqs = length(frequencies);
num_cooling_configs = length(heat_fluxes);
num_decaps = length(decap_ratios);
num_wire_resistivities = length(wire_resistivities);
total_configs = num_stacks * num_perms * num_thicks * num_freqs * num_cooling_configs * num_decaps * num_wire_resistivities;

power = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities);
power_density = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities);

wire_power = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities);
rep_power = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities);
temp = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities);
thickness = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities);
npads = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities);
cap_density = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities);
Ltsv_m2 = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities);


ild_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities);
psn_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities);
power_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities);
wire_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities);
chip_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities);
tsv_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities);

t_sweep_start = cputime;
cur_config = 0;
for cind = 1:num_cooling_configs
    for dind = 1:num_decaps
        for thind = 1:num_thicks
            for nind = 1:length(tiers)
                for pind = 1:num_perms
                    for freq_ind = 1:num_freqs
                        for wire_res_ind = 1:num_wire_resistivities
                            cur_config = cur_config + 1;
                            die_thickness = thicknesses(thind);
                            num_layers_per_block = tiers(nind);
                            epsrd = rel_permittivities(pind);
                            fmax_core = frequencies(freq_ind);
                            wire_resistivity = wire_resistivities(wire_res_ind);

                            fprintf('===============================\n')
                            fprintf('==   cooling: %d/%d \t=====\n',cind,num_cooling_configs);
                            fprintf('==     decap: %d/%d \t=====\n',dind,num_decaps);
                            fprintf('== thickness: %d/%d \t=====\n',thind,num_thicks);
                            fprintf('==     Tiers: %d/%d \t=====\n',nind,num_stacks);
                            fprintf('==     epsrd: %d/%d \t=====\n',pind,num_perms);
                            fprintf('==      freq: %d/%d \t=====\n',freq_ind,num_freqs);
                            fprintf('==      freq: %d/%d \t=====\n',wire_res_ind,num_wire_resistivities);
                            fprintf('==   Overall: %d/%d \t=====\n',cur_config,total_configs);
                            fprintf('===============================\n')

                            %% define parameters

                            [core.chip, core.transistor, core.gate, core.tsv, core.wire, core.psn] = generate_basic_processor_settings(rent_exp_logic,num_layers_per_block,Ng_core,Ach_mm2_core,gate_pitch_core,min_pitch_core,Vdd_core,fmax_core,w_trans);
                            %core.psn.mismatch_tolerance = 0.01;
                            %% Tweak wiring parameters
    %                         core.wire.repeater_fraction = [0.3]; % 1 is default from gen_basic_proc_settings
    %                         core.wire.routing_efficiency = [0.6]; % 0.4 is default from gen_basic_proc_settings
                            core.wire.repeater_fraction = [0.5]; % 1 is default from gen_basic_proc_settings
                            core.wire.routing_efficiency = [0.5]; % 0.4 is default from gen_basic_proc_settings
                            core.wire.repeater_max_area_fraction = 0.3; % (-) Fraction of chip area that can be consumed by repeater/buffer gates
                            core.wire.repeater_via_max_area_fraction = 0.05; % (-) Fraction of routable wire area that can be consumed by vias for repeater connections
                            core.wire.resistivity = wire_resistivity;

                            core.wire.use_graphene = 0;
                            simulation.force_thickness = force_thickness;
                            core.chip.thickness_nominal = die_thickness;
                            core.wire.dielectric_epsr = epsrd;
                            core.psn.decap_area_fraction = decap_ratios(dind);

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

                            heat.up = heat_fluxes(cind);        % above chip

                            %% calculate block parameters
                            [core.chip, core.power, core.tsv, core.wire, core.repeater, core.psn] = codesign_block(core.chip,core.tsv,core.gate,core.transistor,core.wire,heat,core.psn,simulation);

                            power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind) = core.power.total;
                            power_density(cind,dind,thind,nind,pind,freq_ind,wire_res_ind) = core.power.density;

                            wire_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind) = core.power.wiring;
                            rep_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind) = core.power.repeater;
                            temp(cind,dind,thind,nind,pind,freq_ind,wire_res_ind) = core.chip.temperature;
                            thickness(cind,dind,thind,nind,pind,freq_ind,wire_res_ind) = core.chip.thickness;
                            npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind) = core.psn.Npads;


                            ild_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind} = core.chip.iidf;
                            psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind} = core.psn;
                            power_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind} = core.power;
                            wire_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind} = core.wire;
                            chip_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind} = core.chip;
                            tsv_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind} = core.tsv;

                            Ltsv_m2(cind,dind,thind,nind,pind,freq_ind,wire_res_ind) = psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind}.Ltsv/psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind}.l_unit_cell^2;
                            cap_density(cind,dind,thind,nind,pind,freq_ind,wire_res_ind) = psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind}.cap_density;
                        end
                    end
                end
            end
        end
    end
end

t_sweep_stop = cputime;
fprintf('\nTotal time elapsed for parameter sweep: %.3g seconds\n\n',(t_sweep_stop-t_sweep_start));



%% Communication power fraction
% figure(1)
% clf
% hold on
% col = [0 0 0; 0 0 1; 0 1 0 ; 1 0 0];
% for nind = 1:num_stacks
%     pvec = zeros(1,num_freqs);
%     pow_vec = zeros(1,num_freqs);
%     comm_pow_vec = zeros(1,num_freqs);
%     pow_vec(1:end) = power(1,1,1,nind,1,:,1);
%     comm_pow_vec(1:end) = wire_power(1,1,1,nind,1,:,1) + rep_power(1,1,1,nind,1,:,1);
%     comm_pow_frac_vec = comm_pow_vec./pow_vec;
%     plot(frequencies/1e9,comm_pow_frac_vec,'linestyle','-','color',col(nind,:)) ;
% end
% xlim([1 10])
% ylim([0 1])
% xlabel('Clock Frequency (GHz)')
% ylabel('On-chip communication power fraction')
% fixfigs(1,3,14,12)
% %npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind);

%%

figure(1)
clf
hold on
linecol = [ 0 0 0; 0 0 1; 0 1 0; 1 0 0];
for nind = 1:num_stacks
    plvec = zeros(1,num_perms);
    plvec(1:end) = temp(cind,dind,thind,nind,:,freq_ind,wire_res_ind); 
    plot(rel_permittivities,plvec,'color',linecol(nind,:),'linewidth',2);

end

figure(2)
clf
hold on
linecol = [ 0 0 0; 0 0 1; 0 1 0; 1 0 0];
for nind = 1:num_stacks
    plvec = zeros(1,num_perms);
    plvec(1:end) = power(cind,dind,thind,nind,:,freq_ind,wire_res_ind); 
    plot(rel_permittivities,plvec,'color',linecol(nind,:),'linewidth',2);

end


%% Metal pitch vs resistivity
% intel_pitch = [112.5 112.5 112.5 168.8 225 337.6 450.1 566.5 19400  ]; % Actual intel data
% intel_pitch_no_pow = [112.5 112.5 112.5 168.8 225 337.6 450.1 566.5];
% intel_power_total = 73;
% 
% figure(1)
% clf
% hold on
% plot(intel_pitch_no_pow,'k','linewidth',2)
% plot(wire_cell{1,1,1,1,1,6,1}.pn*1e9,'b-','linewidth',2)
% plot(wire_cell{1,1,1,1,1,6,2}.pn*1e9,'r-','linewidth',2)
% xlabel('Wiring tier')
% ylabel('Wire Pitch (nm)')
% 
% figure(2)
% clf
% hold on
% plvec = zeros(1,num_freqs);
%     plvec(1:end) = power(1,1,1,1,1,:,1);
% plot(frequencies/1e9,plvec,'b-','linewidth',2)
% plvec = zeros(1,num_freqs);
%     plvec(1:end) = power(1,1,1,1,1,:,2);
% plot(frequencies/1e9,plvec,'r-','linewidth',2)
% xlabel('Frequency')
% ylabel('Power consumption (W)')
% 
% 
% figure(3)
% clf
% hold on
% plvec = zeros(1,num_freqs);
%     plvec(1:end) = npads(1,1,1,1,1,:,1);
% plot(frequencies/1e9,plvec,'b-','linewidth',2)
% plvec = zeros(1,num_freqs);
%     plvec(1:end) = npads(1,1,1,1,1,:,2);
% plot(frequencies/1e9,plvec,'r-','linewidth',2)
% xlabel('Frequency')
% ylabel('Number of power TSVs')


