%% Simulation parameters
simulation.skip_psn_loops = 1; % Skip PSN TSV homing for faster debug
simulation.skip_thermal = 1; % Skip thermal analysis for faster debug

simulation.use_joyner = 0;
simulation.redo_wiring_after_repeaters = 0;
simulation.topdown_WLARI = 1; % Use topdown simultaneous WLA and RI (0 = use standard bottom-up optimal WLA, followed by one pass of RI)
simulation.separate_wiring_tiers = 1; % 1 = Each logic plane will have its own wiring tiers between it and the next logic plane
                                      % 0 = All metal layers for entire device will be routed on top of entire 3D stack

simulation.draw_thermal_map = 0; % Plot thermal profile of each chip
simulation.print_thermal_data = 0; % Output max temp in each layer to console

simulation.wla_max_attempts = 30; % 15 is default
simulation.wla_min_bot_fill_factor = 0.90; % 0.97 is default
simulation.wla_min_top_fill_factor = 0.01; % 0.01 is default

simulation.freq_binsearch = 0;
simulation.freq_binsearch_target = 90;
simulation.freq_binsearch_raw_tol = 0.05;
simulation.freq_binsearch_max_gens = 10;


%% Logic core parameters
compression_factor = 1; % linear scaling factor. 1 = actual 32nm design, 4.57 = equivalent 7nm SB
Ng_core = 86e6/4; %86M transistors, assume 2in NAND -> /4 to get total NAND gates
Ach_mm2_core = 18.5/(compression_factor^2);
gate_pitch_core = 465e-9*2/compression_factor;
min_pitch_core = 112.5e-9/compression_factor;
fmax_core = 3.5e9;
w_trans = 32e-9/compression_factor;
Vdd_core = 1.25;

%% Thermal parameters

r_air = 1/1.825; %K/W for a 1cm^2 HS % alt, 0.6
r_water = 1/4.63; %K/W for a 1cm^2 HS
A_hs = (1e-2)^2; % 1 cm^2

h_air = 1/(r_air*A_hs);
h_water = 1/(r_water*A_hs);


%%
rent_exp_logic = 0.6;
rent_exp_mem = 0.4;
rent_exp_gpu = 0.55;


%%
rho_ag = 15.9e-9;
rho_au = 24.4e-9;
rho_cu = 17.2e-9;
rho_w = 56.0e-9;
rho_ni = 69.9e-9;
%rho_co_al = 0e-9;
rho_al = 26.5e-9;
rho_all_mets = [rho_ag rho_cu rho_au rho_al rho_w rho_ni];


%% 
tiers = [1 2 4 8];
thicknesses = logspace(-6,-4,51);
force_thickness = 1;
rel_permittivities = [3.0];
frequencies = fmax_core; % if simulation.freq_binsearch is set, (1) is min freq and (2) is max freq
heat_fluxes = [ h_air];
decap_ratios = [0.1];%[0.01 0.1 1];
%wire_resistivities = [rho_ag rho_cu rho_au rho_al rho_w rho_ni];
wire_resistivities = [rho_cu];
wire_material_flags = {'00'}; % binary strings. bit1 = use_graphene, bit0 = use alt_em_mat
scaling_factor = [1];


num_stacks = length(tiers);
num_perms = length(rel_permittivities);
num_thicks = length(thicknesses);

if (simulation.freq_binsearch == 1)
    num_freqs = 1; % skip freq loops
else
    num_freqs = length(frequencies);
end
num_cooling_configs = length(heat_fluxes);
num_decaps = length(decap_ratios);
num_wire_resistivities = length(wire_resistivities);
num_wire_flags = length(wire_material_flags);
num_scaling_factors = length(scaling_factor);
total_configs = num_stacks * num_perms * num_thicks * num_freqs * num_cooling_configs * num_decaps * num_wire_resistivities * num_wire_flags * num_scaling_factors;

power = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
power_density = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
freq = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
wire_power = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
rep_power = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
temp = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
thickness = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
npads = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
cap_density = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
Ltsv_m2 = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);


ild_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
psn_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
power_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
wire_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
chip_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
tsv_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);

t_sweep_start = cputime;
cur_config = 0;
for cind = 1:num_cooling_configs
    for dind = 1:num_decaps
        for thind = 1:num_thicks
            for nind = 1:length(tiers)
                for pind = 1:num_perms
                    for freq_ind = 1:num_freqs
                        for wire_res_ind = 1:num_wire_resistivities
                            for wire_flag_ind = 1:num_wire_flags
                                for scaling_ind = 1:num_scaling_factors
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
                                    fprintf('==      freq: %d/%d \t=====\n',freq_ind,num_freqs);
                                    fprintf('==       rho: %d/%d \t=====\n',wire_res_ind,num_wire_resistivities);
                                    fprintf('==       mat: %d/%d \t=====\n',wire_flag_ind,num_wire_flags);
                                    fprintf('==     scale: %d/%d \t=====\n',scaling_ind,num_scaling_factors);
                                    fprintf('==   Overall: %d/%d \t=====\n',cur_config,total_configs);
                                    fprintf('===============================\n')
                                    
                                    %% Rescale dimensions if necessary
                                    compression_factor = scaling_factor(scaling_ind); % linear scaling factor. 1 = actual 32nm design, 4.57 = equivalent 7nm SB
                                    Ach_mm2_scaled = Ach_mm2_core/compression_factor^2;
                                    gate_pitch_scaled = gate_pitch_core/compression_factor;
                                    min_pitch_scaled = min_pitch_core/compression_factor;
                                    w_trans_scaled = w_trans/compression_factor;
                                    

                                    %% define parameters
                                    [core.chip, core.transistor, core.gate, core.tsv, core.wire, core.psn, core.heat] = generate_basic_processor_settings(rent_exp_logic,num_layers_per_block,Ng_core,Ach_mm2_scaled,gate_pitch_scaled,min_pitch_scaled,Vdd_core,fmax_core,w_trans_scaled);
                                    %core.psn.mismatch_tolerance = 0.01;
                                    %% Tweak wiring parameters
            %                         core.wire.repeater_fraction = [0.3]; % 1 is default from gen_basic_proc_settings
            %                         core.wire.routing_efficiency = [0.6]; % 0.4 is default from gen_basic_proc_settings
                                    core.wire.repeater_fraction = [0.4]; % 1 is default from gen_basic_proc_settings
                                    core.wire.routing_efficiency = [0.5]; % 0.4 is default from gen_basic_proc_settings
                                    core.wire.repeater_max_area_fraction = 0.2; % (-) Fraction of chip area that can be consumed by repeater/buffer gates
                                    core.wire.repeater_via_max_area_fraction = 0.05; % (-) Fraction of routable wire area that can be consumed by vias for repeater connections
                                    core.gate.output_resistance = 8e3*compression_factor; % Ohm
                                    core.transistor.capacitance = 1e-15*1e6*3*w_trans; % ITRS projection is 1fF/um of gate width. This is an estimate for pMOS transistor capacitance
                                    core.transistor.leakage_current_per_micron = 100e-9; %(A/um) % 32nm IOFF

                                    wire_flag_key = bin2dec(wire_material_flags(wire_flag_ind));
                                    alt_met_flag = bitget(wire_flag_key,1);
                                    graphene_flag = bitget(wire_flag_key,2);
                                    core.wire.use_graphene = graphene_flag;
                                    core.wire.use_em_resistant_metal = alt_met_flag;   % (-) Allow or disallow use of electromigration-resistant metals below a specified minimum pitch
                                    core.wire.min_non_em_width = 25e-9; % (m) If use_em_resistant_metal is set to 1, Cu resistivity will be replaced with wire.alt_resistivity_em below this pitch
                                    if(alt_met_flag == 1)
                                        core.wire.alt_resistivity_em = wire_resistivity;
                                        core.wire.resistivity = rho_cu;
                                    else
                                        core.wire.resistivity = wire_resistivity;
                                    end
                                    
                                    simulation.force_thickness = force_thickness;
                                    core.chip.thickness_nominal = die_thickness;
                                    core.wire.dielectric_epsr = epsrd;
                                    core.psn.decap_area_fraction = decap_ratios(dind);

                                    core.heat.up = heat_fluxes(cind);        % above chip
                                    
                                    %% If we're searching for frequencies below a certain temperature
                                    if (simulation.freq_binsearch == 1)
                                        fmin = frequencies(1);
                                        fmax = frequencies(2);
                                        max_gens = simulation.freq_binsearch_max_gens;
                                        target_max_value = simulation.freq_binsearch_target;
                                        target_cur_value = target_max_value*2; % start with something invalid so we run it at least once
                                        abs_err = abs(target_max_value - target_cur_value);
                                        tolerance = simulation.freq_binsearch_raw_tol;
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
                                    else
                                        %% calculate block parameters
                                        [core.chip, core.power, core.tsv, core.wire, core.repeater, core.psn] = codesign_block(core.chip,core.tsv,core.gate,core.transistor,core.wire,core.heat,core.psn,simulation);
                                    end
                                    power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.power.total;
                                    power_density(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.power.density;
                                    freq(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = 1/core.chip.clock_period;

                                    wire_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.power.wiring;
                                    rep_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.power.repeater;
                                    temp(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.chip.temperature;
                                    thickness(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.chip.thickness;
                                    npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.psn.Npads;


                                    ild_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind} = core.chip.iidf;
                                    psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind} = core.psn;
                                    power_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind} = core.power;
                                    wire_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind} = core.wire;
                                    chip_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind} = core.chip;
                                    tsv_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind} = core.tsv;

                                    if (simulation.skip_psn_loops == 0)
                                        Ltsv_m2(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind}.Ltsv/psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind}.l_unit_cell^2;
                                        cap_density(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind}.cap_density;
                                    end
                                end
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

figure(4)
clf
hold on
linecol = [ 0 0 0; 0 0 1; 0 1 0; 1 0 0];
for nind = 1:num_stacks
    plvec = zeros(1,num_thicks);
    plvec(1:end) = wire_power(cind,dind,:,nind,pind,freq_ind,wire_res_ind) + rep_power(cind,dind,:,nind,pind,freq_ind,wire_res_ind);
    plot(thicknesses/1e-6,plvec,'linestyle','-','color',linecol(nind,:),'linewidth',2);
end
set(gca,'xscale','log')
h = xlabel('Tier thickness (microns)');
ylabel('On-chip comm. power (W)')
fixfigs(4,2,14,12)

%% Communication power fraction
figure(1)
clf
hold on
col = [0 0 0; 0 0 1; 0 1 0 ; 1 0 0];
for nind = 1:num_stacks
    pvec = zeros(1,num_freqs);
    pow_vec = zeros(1,num_freqs);
    comm_pow_vec = zeros(1,num_freqs);
    pow_vec(1:end) = power(1,1,1,nind,1,:,1);
    comm_pow_vec(1:end) = wire_power(1,1,1,nind,1,:,1) + rep_power(1,1,1,nind,1,:,1);
    comm_pow_frac_vec = comm_pow_vec./pow_vec;
    plot(frequencies/1e9,comm_pow_frac_vec,'linestyle','-','color',col(nind,:)) ;
end
xlim([1 10])
ylim([0 1])
xlabel('Clock Frequency (GHz)')
ylabel('On-chip communication power fraction')
fixfigs(1,2,14,12)
%npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_mat_ind,scaling_ind);

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
% nind = [1 2 4 8];
% area = zeros(1,num_stacks);
% area_mm2 = zeros(1,num_stacks);
% for nnn = 1:4
%     plvec = zeros(1,num_decaps);
%     plvec(1:end) = npads(cind,:,thind,(nnn),pind,freq_ind,wire_res_ind);
%     area(nnn) = chip_cell{cind,dind,thind,nnn,pind,freq_ind,wire_res_ind}.area_per_layer_m2;
%     area_mm2(nnn) = area(nnn)/1e-6;
%     plot(decap_ratios,plvec,'color',linecol(nnn,:),'linewidth',2);
% end
% set(gca,'yscale','log')
% set(gca,'xscale','log')
% xlabel('Fraction of die used for decoupling capacitors')
% ylabel('Number of power and ground TSVs')
% fixfigs(2,2,14,12)
% 
% figure(3)
% clf
% hold on
% linecol = [ 0 0 0; 0 0 1; 0 1 0; 1 0 0];
% nind = [1 2 4 8];
% area = zeros(1,num_stacks);
% area_mm2 = zeros(1,num_stacks);
% for nnn = 1:4
%     plvec = zeros(1,num_decaps);
%     plvec(1:end) = npads(cind,:,thind,(nnn),pind,freq_ind,wire_res_ind);
%     area(nnn) = chip_cell{cind,dind,thind,nnn,pind,freq_ind,wire_res_ind}.area_per_layer_m2;
%     area_mm2(nnn) = area(nnn)/1e-6;
%     plot(decap_ratios*area_mm2(nnn),plvec,'color',linecol(nnn,:),'linewidth',2);
% end
% set(gca,'yscale','log')
% set(gca,'xscale','log')
% xlabel('Area used for decoupling capacitors (mm^2)')
% ylabel('Number of power and ground TSVs')
% fixfigs(3,2,14,12)

%% Power TSVs, Tiers, and ILD

% Use these parameters to generate this plot
% tiers = [1 2 4 8];
% thicknesses = [1e-6 10e-6 100e-6];
% force_thickness = 1;
% rel_permittivities = [1 3.9];
% frequencies = fmax_core;
% heat_fluxes = [ h_air ];
% decap_ratios = [0.1];%[0.01 0.1 1];
% wire_resistivities = [17.2e-9];

% num_ptsvs = zeros(num_stacks,2*num_perms);
% for nind = 1:num_stacks
%     for pind = 1:num_perms
%         th = 0;
%         for thind = [1 3]
%             th = th + 1;
%             num_ptsvs(nind,(pind-1)*2 + th) = npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_mat_ind,scaling_ind);           
%         end
%     end
% end
% 
% figure(1)
% clf
% b = bar(num_ptsvs,1,'grouped');
% colormap jet
% set(gca,'xticklabel',{'1','2','4','8'})
% ylim([1e0 1e4])
% set(gca,'yscale','log')
% xlabel('Number of tiers')
% ylabel('Number of power delivery TSVs')
% b(1).FaceColor = 'blue';
% b(2).FaceColor = 'green';
% b(3).FaceColor = 'yellow';
% b(4).FaceColor = 'red';
% fixfigs(1,2,14,12)

%% GNRs in scaled and unscaled core
% 
% figure(1)
% clf
% hold all
% for scaling_ind = 1:num_scaling_factors
%     for wire_mat_ind = 1:num_wire_flags
%         wire_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_mat_ind,scaling_ind};
%         plot(wire_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_mat_ind,scaling_ind}.pn*1e9);
%     end
% end
% fixfigs(1,2,14,12)

%% Impact of wire resistivity
% figure(1)
% clf
% hold all
% pvec = zeros(1,num_wire_resistivities);
% for wire_res_ind = 1:num_wire_resistivities
%     pvec(wire_res_ind) = length(wire_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,1,1}.pn);
% end
% plot(wire_resistivities*1e9,pvec,'k-');
% 
% pvec = zeros(1,num_wire_resistivities);
% for wire_res_ind = 1:num_wire_resistivities
%     pvec(wire_res_ind) = length(wire_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,1,2}.pn);
% end
% plot(wire_resistivities*1e9,pvec,'b-');
% 
% pvec = zeros(1,num_wire_resistivities);
% for wire_res_ind = 1:num_wire_resistivities
%     pvec(wire_res_ind) = length(wire_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,2,2}.pn);
% end
% plot(wire_resistivities*1e9,pvec,'r--');
% fixfigs(1,2,14,12)
% 
% wire_res_ind = 2;
% intel_pitch_no_pow = [112.5 112.5 112.5 168.8 225 337.6 450.1 566.5];
% figure(2)
% clf
% hold on
% plot(intel_pitch_no_pow,'k')
% plot(wire_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,1,1}.pn*1e9)
% fixfigs(2,2,14,12)


%% Power efficiency vs frequency
% % Inputs
% % tiers = [1 2 4 8];
% % thicknesses = [10e-6];
% % force_thickness = 1;
% % rel_permittivities = 3.0;
% % frequencies = linspace(0.1e9, 5e9, 1e2);
% % heat_fluxes = [ h_air ];
% % decap_ratios = [0.1];%[0.01 0.1 1];
% % %wire_resistivities = [rho_ag rho_cu rho_au rho_al rho_w rho_ni];
% % wire_resistivities = [rho_cu];
% % wire_material_flags = {'00'}; % binary strings. bit1 = use_graphene, bit0 = use alt_em_mat
% % scaling_factor = [1];
% 
% % Plots
% figure(1)
% clf
% hold on
% for nind = 1:num_stacks
%     colors = [ 0 0 0 ; 0 0 1; 0 1 0 ; 1 0 0 ];
%     pow_tot = zeros(1,num_freqs);
%     pow_wire = zeros(1,num_freqs);
%     pow_rep = zeros(1,num_freqs);
%     
%     pow_tot(1,:) = power(cind,dind,thind,nind,pind,:,wire_res_ind,wire_flag_ind,scaling_ind);
%     pow_wire(1,:) = wire_power(cind,dind,thind,nind,pind,:,wire_res_ind,wire_flag_ind,scaling_ind);
%     pow_rep(1,:) = rep_power(cind,dind,thind,nind,pind,:,wire_res_ind,wire_flag_ind,scaling_ind);
%     pow_comm = pow_wire + pow_rep;
%     
%     
%     pow_logic = pow_tot - (pow_wire + pow_rep);
%     pow_eff = pow_logic./pow_tot;
%     pow_comm_ratio = pow_comm./pow_tot;
%     pow_log_ratio = 1 - pow_comm_ratio;
%     
%     plot(frequencies/1e9,pow_comm_ratio,'color',colors(nind,:))
%     xlim([0.5 5])
%     xlabel('Clock Frequency (GHz)')
%     ylabel('On-chip Communication Power Fraction')
%     %set(gca,'xscale','log')
%     fixfigs(1,2,14,12)
% end
%     
%     
% %     
% % power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.power.total;
% % power_density(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.power.density;
% % 
% % wire_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.power.wiring;
% % rep_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.power.repeater;


%%  Max frequency vs ILD
% % Inputs
% tiers = [1 2 4 8];
% thicknesses = [10e-6];
% force_thickness = 1;
% rel_permittivities = linspace(1,4,41);
% frequencies = [3e8 1e10]; % if simulation.freq_binsearch is set, (1) is min freq and (2) is max freq
% heat_fluxes = [ h_air h_water];
% decap_ratios = [0.1];%[0.01 0.1 1];
% %wire_resistivities = [rho_ag rho_cu rho_au rho_al rho_w rho_ni];
% wire_resistivities = [rho_cu];
% wire_material_flags = {'00'}; % binary strings. bit1 = use_graphene, bit0 = use alt_em_mat
% scaling_factor = [1];


% Plots
% figure(1)
% clf
% hold on
% 
% for nind = 1:num_stacks
%     colors = [ 0 0 0 ; 0 0 1; 0 1 0 ; 1 0 0 ];
%     fr_vec1 = zeros(1,num_perms);
%     fr_vec2 = zeros(1,num_perms);
%     cind = 1;
%     fr_vec1(1,:) = freq(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind)/1e9;
%     cind = 2;
%     fr_vec2(1,:) = freq(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind)/1e9;
%     
%     plot(rel_permittivities,fr_vec1,'color',colors(nind,:),'linestyle','-')
%     plot(rel_permittivities,fr_vec2,'color',colors(nind,:),'linestyle','--')
% end
% xlabel('ILD Relative Permittivity')
% ylabel('Maximum Frequency')
% %set(gca,'yscale','log')
% fixfigs(1,2,14,12)
% 
% figure(2)
% clf
% hold on
% p_ave_vec = zeros(1,num_stacks);
% p_ave_vec2 = zeros(1,num_stacks);
% comm_pow_ave_vec = zeros(1,num_stacks);
% comm_pow_ave_vec2 = zeros(1,num_stacks);
% %pfrac_vec = zeros(1,num_stacks);
% %pfrac_vec2 = zeros(1,num_stacks);
% for nind = 1:num_stacks
%     colors = [ 0 0 0 ; 0 0 1; 0 1 0 ; 1 0 0 ];
%     cind = 1;
%     pow_vec = zeros(1,num_perms);
%     comm_pow_vec = zeros(1,num_perms);
%     pow_vec(1,:) = power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind);
%     comm_pow_vec(1:end) = wire_power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) + rep_power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind);
%     p_ave_vec(nind) = mean(pow_vec);
%     comm_pow_ave_vec(nind) = mean(comm_pow_vec);
%     
%     
%     cind = 2;
%     pow_vec2 = zeros(1,num_perms);
%     comm_pow_vec2 = zeros(1,num_perms);
%     pow_vec2(1,:) = power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind);
%     comm_pow_vec2(1:end) = wire_power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) + rep_power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind);
%     p_ave_vec2(nind) = mean(pow_vec2);
%     comm_pow_ave_vec2(nind) = mean(comm_pow_vec2);
%     
%     
%     %plot(rel_permittivities,comm_pow_vec./pow_vec,'color',colors(nind,:))
%     plot(rel_permittivities,pow_vec,'color',colors(nind,:),'linestyle','-')
%     plot(rel_permittivities,pow_vec2,'color',colors(nind,:),'linestyle','--')
% end
% xlabel('ILD Relative Permittivity')
% ylabel('Power (W)')
% %set(gca,'xscale','log')
% fixfigs(2,2,14,12)
% 
% 
% pfrac_vec = comm_pow_ave_vec./p_ave_vec;
% pfrac_vec2 = comm_pow_ave_vec2./p_ave_vec2;
% 
% figure(3)
% clf
% hold on
% plot(tiers,p_ave_vec,'r')
% plot(tiers,p_ave_vec2,'b')
% xlabel('Tiers')
% ylabel('Power (W)')
% fixfigs(3,2,14,12)
% 
% 
% pmat = [p_ave_vec ; p_ave_vec2]';
% figure(4)
% b = bar(pmat,1,'grouped');
% colormap jet
% set(gca,'xticklabel',{'1','2','4','8'})
% %ylim([1e0 1e4])
% %set(gca,'yscale','log')
% xlabel('Number of tiers')
% ylabel('Power (W)')
% b(1).FaceColor = 'blue';
% % b(2).FaceColor = 'green';
% b(2).FaceColor = 'yellow';
% %b(2).FaceColor = 'red';
% fixfigs(4,2,14,12)
% 
% pmat = [pfrac_vec ; pfrac_vec2]';
% figure(5)
% b = bar(pmat,1,'grouped');
% colormap jet
% set(gca,'xticklabel',{'1','2','4','8'})
% %ylim([1e0 1e4])
% %set(gca,'yscale','log')
% xlabel('Number of tiers')
% ylabel('Comm Power Fraction')
% b(1).FaceColor = 'blue';
% % b(2).FaceColor = 'green';
% b(2).FaceColor = 'yellow';
% %b(2).FaceColor = 'red';
% fixfigs(5,2,14,12)
% 
% 
% % Plots
% figure(6)
% clf
% hold on
% inds = [1 14 28 41];
% for nind = 1:num_stacks
%     colors = [ 0 0 0 ; 0 0 1; 0 1 0 ; 1 0 0 ];
%     fr_vec1 = zeros(1,num_perms);
%     fr_vec2 = zeros(1,num_perms);
%     pow_vec = zeros(1,num_perms);
%     pow_vec2 = zeros(1,num_perms);
%     cind = 1;
%     fr_vec1(1,:) = freq(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind);
%     pow_vec(1,:) = power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind);
%     epc_vec = pow_vec./fr_vec1;
%     
%     cind = 2;
%     fr_vec2(1,:) = freq(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind);
%     pow_vec2(1,:) = power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind);
%     epc_vec2 = pow_vec2./fr_vec2;
%     
%     
%     plot(rel_permittivities,epc_vec*1e9,'color',colors(nind,:),'linestyle','-')
%     plot(rel_permittivities,epc_vec2*1e9,'color',colors(nind,:),'linestyle','--')
% end
% xlabel('ILD Relative Permittivity')
% ylabel('Energy per cycle (nJ)')
% %set(gca,'yscale','log')
% fixfigs(6,2,14,12)
% 
% 
% epc_mat = zeros(4,num_stacks);

%% Power consumption vs power density

% % Inputs
% tiers = 1:8;
% thicknesses = [10e-6];
% force_thickness = 1;
% rel_permittivities = [3.0];
% frequencies = fmax_core; % if simulation.freq_binsearch is set, (1) is min freq and (2) is max freq
% heat_fluxes = [ h_air];
% decap_ratios = [0.1];%[0.01 0.1 1];
% %wire_resistivities = [rho_ag rho_cu rho_au rho_al rho_w rho_ni];
% wire_resistivities = [rho_cu];
% wire_material_flags = {'00'}; % binary strings. bit1 = use_graphene, bit0 = use alt_em_mat
% scaling_factor = [1];


pvec = zeros(1,num_stacks);
pdens_vec = zeros(1,num_stacks);
pvec(1,:) = power(cind,dind,thind,:,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind);
pdens_vec(1,:) = power_density(cind,dind,thind,:,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind)/100^2;

figure(1)
clf
[ax, h1, h2] = plotyy(tiers,pvec,tiers,pdens_vec);
set(ax(1),'ycolor','k')
set(ax(2),'ycolor','k')
ax(1).FontSize = 12;
ax(2).FontSize = 12;
%ax(1).YLabel.String = 'Power (W)';
%ax(2).YLabel.String = 'Power Density (W/cm^2)';
xlabel('Number of tiers')
fixfigs(1,2,14,12)



