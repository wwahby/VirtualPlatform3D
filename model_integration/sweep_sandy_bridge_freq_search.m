%% Simulation parameters
simulation.use_joyner = 0;
simulation.redo_wiring_after_repeaters = 0;
simulation.topdown_WLARI = 1; % Use topdown simultaneous WLA and RI (0 = use standard bottom-up optimal WLA, followed by one pass of RI)
simulation.skip_psn_loops = 1; % Skip PSN TSV homing for faster debug
simulation.draw_thermal_map = 0; % Plot thermal profile of each chip
simulation.print_thermal_data = 0; % Output max temp in each layer to console
simulation.separate_wiring_tiers = 1; % 1 = Each logic plane will have its own wiring tiers between it and the next logic plane
                                      % 0 = All metal layers for entire device will be routed on top of entire 3D stack
simulation.wla_max_attempts = 15;
simulation.wla_min_top_fill_factor = 0.01;

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

r_air = 1/1.825; %K/W for a 1cm^2 HS
r_water = 1/4.63; %K/W for a 1cm^2 HS
A_hs = (1e-2)^2; % 1 cm^2

h_air = 1/(r_air*A_hs);
h_water = 1/(r_water*A_hs);


%%
rent_exp_logic = 0.6;
rent_exp_mem = 0.4;
rent_exp_gpu = 0.55;

%% 
tiers = 1:8;
thicknesses = [10e-6];
force_thickness = 1;
rel_permittivities = linspace(1,4,41);
frequencies = [3e8 1e10]; % (1) is min freq (2) is max freq
heat_fluxes = [ h_air h_water];
decap_ratios = [0.1];%[0.01 0.1 1];
wire_resistivities = [17.2e-9];


num_stacks = length(tiers);
num_perms = length(rel_permittivities);
num_thicks = length(thicknesses);
num_cooling_configs = length(heat_fluxes);
num_decaps = length(decap_ratios);
num_wire_resistivities = length(wire_resistivities);
num_freqs = 1;
total_configs = num_stacks * num_perms * num_thicks * num_cooling_configs * num_freqs * num_decaps * num_wire_resistivities;

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

                            fprintf('\n===============================\n')
                            fprintf('==   cooling: %d/%d \t=====\n',cind,num_cooling_configs);
                            fprintf('==     decap: %d/%d \t=====\n',dind,num_decaps);
                            fprintf('== thickness: %d/%d \t=====\n',thind,num_thicks);
                            fprintf('==     Tiers: %d/%d \t=====\n',nind,num_stacks);
                            fprintf('==     epsrd: %d/%d \t=====\n',pind,num_perms);
                            fprintf('==       rho: %d/%d \t=====\n',wire_res_ind,num_wire_resistivities);
                            fprintf('==   Overall: %d/%d \t=====\n',cur_config,total_configs);
                            fprintf('===============================\n')

                            %% define parameters

                            [core.chip, core.transistor, core.gate, core.tsv, core.wire, core.psn, core.heat] = generate_basic_processor_settings(rent_exp_logic,num_layers_per_block,Ng_core,Ach_mm2_core,gate_pitch_core,min_pitch_core,Vdd_core,fmax_core,w_trans);
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

                            core.heat.up = heat_fluxes(cind);        % above chip

                            %% calculate block parameters
                            fmin = frequencies(1);
                            fmax = frequencies(2);
                            max_gens = 10;
                            target_max_value = 90;
                            target_cur_value = target_max_value*2; % start with something invalid so we run it at least once
                            abs_err = abs(target_max_value - target_cur_value);
                            tolerance = 0.5;
                            left = fmin;
                            right = fmax;
                            
                            num_gens = 0;
                            time_bin_start = cputime;
                            while ((abs_err > tolerance) && (num_gens < max_gens))
                                mid = 1/2*(left+right);
                                core.chip.clock_period = 1/mid;
                                [core.chip, core.power, core.tsv, core.wire, core.repeater, core.psn] = codesign_block(core.chip,core.tsv,core.gate,core.transistor,core.wire,core.heat,core.psn,simulation);
                                target_cur_value = core.chip.temperature;
                                abs_err = abs(target_max_value - target_cur_value);

                                if (target_cur_value > target_max_value)
                                    right = mid;
                                elseif (target_cur_value < target_max_value)
                                    left = mid;
                                else
                                    left = 1/2*(left + mid);
                                    right = 1/2*(right + mid);
                                end

                                num_gens = num_gens + 1;
                                fprintf('Run %d: \t Freq: %.3g \t Temp: %.4d\n\n',num_gens, mid, target_cur_value);
                            end
                            time_bin_stop = cputime;
                            fprintf('Time for last binary search: %d seconds\n',time_bin_stop - time_bin_start);

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

                            if (simulation.skip_psn_loops == 0)
                                Ltsv_m2(cind,dind,thind,nind,pind,freq_ind,wire_res_ind) = psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind}.Ltsv/psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind}.l_unit_cell^2;
                                cap_density(cind,dind,thind,nind,pind,freq_ind,wire_res_ind) = psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind}.cap_density;
                            end
                        end
                    end
                end
            end
        end
    end
end

t_sweep_stop = cputime;
fprintf('\nTotal time elapsed for parameter sweep: %.3g seconds\t(%.3g minutes)\n\n',(t_sweep_stop-t_sweep_start),(t_sweep_stop-t_sweep_start)/60);


%% Power vs tier thickness (tier thickness sweep only)
% 
% figure(1)
% clf
% hold on
% linecol = [ 0 0 0; 0 0 1; 0 1 0; 1 0 0];
% for nind = 1:num_stacks
%     plvec = zeros(1,num_thicks);
%     plvec(1:end) = wire_power(cind,dind,:,nind,pind,freq_ind,wire_res_ind) + rep_power(cind,dind,:,nind,pind,freq_ind,wire_res_ind);
%     plot(thicknesses/1e-6,plvec,'linestyle','-','color',linecol(nind,:),'linewidth',2);
% end
% set(gca,'xscale','log')
% h = xlabel('Tier thickness (microns)');
% ylabel('On-chip communication power (W)')
% fixfigs(1,2,14,12)

% %% Communication power fraction
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
% fixfigs(1,2,14,12)
% %npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind);

%% Use these for tier and ILD sweep
% 
% figure(1)
% clf
% hold on
% linecol = [ 0 0 0; 0 0 1; 0 1 0; 1 0 0];
% nn = 0;
% for nind = [1 2 4 8]
%     nn = nn+1;
%     plvec = zeros(1,num_perms);
%     plvec(1:end) = temp(cind,dind,thind,nind,:,freq_ind,wire_res_ind); 
%     plot(rel_permittivities,plvec,'color',linecol(nn,:),'linewidth',2);
% 
% end
% xlabel('ILD Relative Permittivity')
% ylabel('Maximum Temperature (C)')
% 
% figure(2)
% clf
% hold on
% linecol = [ 0 0 0; 0 0 1; 0 1 0; 1 0 0];
% nn = 0;
% for nind = [1 2 4 8]
%     nn = nn + 1;
%     plvec = zeros(1,num_perms);
%     plvec(1:end) = power(cind,dind,thind,nind,:,freq_ind,wire_res_ind); 
%     plot(rel_permittivities,plvec,'color',linecol(nn,:),'linewidth',2);
% end
% xlabel('ILD Relative Permittivity')
% ylabel('Power Consumption (W)')
% 
% figure(7)
% clf
% hold on
% linecol = [ 0 0 0; 0 0 1; 0 1 0; 1 0 0];
% for pind = 1:num_perms
%     plvec = zeros(1,num_stacks);
%     plvec(1:end) = power(cind,dind,thind,:,pind,freq_ind,wire_res_ind); 
%     plot(tiers,plvec,'color',linecol(5-pind,:),'linewidth',2);
% end
% xlabel('Number of tiers')
% ylabel('Power Consumption (W)')
% 
% 
% figure(3)
% clf
% hold on
% plvec = zeros(1,num_perms);
% avec = zeros(3,num_perms);
% plvec(1:end) = power(cind,dind,thind,nind,:,freq_ind,wire_res_ind);
% avec(1,:) = plvec;
% plot(rel_permittivities,plvec,'k','linewidth',2);
% plvec(1:end) = wire_power(cind,dind,thind,nind,:,freq_ind,wire_res_ind);
% avec(2,:) = plvec;
% plot(rel_permittivities,plvec,'b','linewidth',2);
% plvec(1:end) = rep_power(cind,dind,thind,nind,:,freq_ind,wire_res_ind);
% avec(3,:) = plvec;
% avec(1,:) = avec(1,:) - avec(2,:) - avec(3,:);
% plot(rel_permittivities,plvec,'r','linewidth',2);
% xlabel('ILD Relative Permittivity')
% ylabel('Single Core Power Consumption (W)')
% fixfigs(3,2,14,12)
% 
% 
% bvec = avec(2:end,:);
% 
% figure(4)
% clf
% hold on
% h = area(rel_permittivities,avec');
% h(1).FaceColor = [0.3 0.3 1];
% h(2).FaceColor = [1 1 0.3];
% h(2).FaceColor = [0.3 1 0.3];
% h(3).FaceColor = [1 0.3 0.3];
% xlabel('ILD Relative Permittivity','fontsize',14)
% ylabel('Single Core Power Consumption (W)','fontsize',14)
% 
% figure(5)
% clf
% hold on
% h = area(rel_permittivities,bvec');
% h(1).FaceColor = [0.3 0.3 1];
% h(2).FaceColor = [1 0.3 0.3];
% xlabel('ILD Relative Permittivity','fontsize',14)
% ylabel('Single Core Power Consumption (W)','fontsize',14)


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


%% Impact of decap

% figure(1)
% clf
% hold on
% linecol = [ 0 0 0; 0 0 1; 0 1 0; 1 0 0];
% for ddd = 1:num_decaps
%     plvec = zeros(1,num_stacks);
%     plvec(1:end) = npads(cind,ddd,thind,:,pind,freq_ind,wire_res_ind);
%     plot(tiers,plvec,'color',linecol(ddd,:),'linewidth',2);
% end
% %set(gca,'yscale','log')
% xlabel('Number of tiers')
% ylabel('Power and ground connections per tier')
% fixfigs(1,2,14,12)
% 
% figure(2)
% clf
% hold on
% linecol = [ 0 0 0; 0 0 1; 0 1 0; 1 0 0];
% for nind= 1:num_stacks
%     plvec = zeros(1,num_decaps);
%     plvec(1:end) = npads(cind,:,thind,nind,pind,freq_ind,wire_res_ind);
%     plot(decap_ratios,plvec,'color',linecol(nind,:),'linewidth',2);
% end
% %set(gca,'yscale','log')
% set(gca,'xscale','log')
% xlabel('Decap ratio')
% ylabel('Number of power and ground TSVs')
% fixfigs(2,2,14,12)

%% Max frequency under 90C

max_freqs = zeros(num_stacks,num_perms);
temp_vec = zeros(num_stacks,num_perms);
cind = 1;
for nind = 1:num_stacks
    for pind = 1:num_perms
        clock_period = chip_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind}.clock_period;
        max_freqs(nind,pind) = 1/clock_period;
        temp_vec(nind,pind) = temp(cind,dind,thind,nind,pind,freq_ind,wire_res_ind);
        pg_vec(nind,pind) = npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind);
    end
end