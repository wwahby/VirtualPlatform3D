%function sweep_data = sweep_design(design, sweep, simulation)
% Just use as script for now
t_sweep_start = cputime;
%% Unpack design settings
rent_exp = design.rent_exp;
Ng_core = design.Ng_core;
Ach_mm2 = design.Ach_mm2;
gate_pitch = design.gate_pitch;
min_pitch = design.min_pitch;
w_trans = design.w_trans;

% %% Unpack sweep settings
% tiers = tiers;
% thicknesses = thicknesses;
% force_thickness = force_thickness;
% rel_permittivities = rel_permittivities;
% frequencies = frequencies;
% heat_fluxes = heat_fluxes;
% decap_ratios = decap_ratios;
% wire_resistivities = wire_resistivities;
% wire_material_flags = wire_material_flags;
% scaling_factors = scaling_factors;
% barrier_thicknesses = barrier_thicknesses;
% barrier_resistivities = barrier_resistivities;

%% Initialize sweep variables
num_stacks = length(tiers);
num_perms = length(rel_permittivities);
num_thicks = length(thicknesses);

if (simulation.freq_binsearch == 1)
    num_freqs = 1; % skip freq loops
else
    num_freqs = length(frequencies);
end

if (simulation.force_power == 1)
    num_forced_powers = length(power_forced_vec);
else
    num_forced_powers = 1;
end

num_cooling_configs = length(heat_fluxes);
num_decaps = length(decap_ratios);
num_wire_resistivities = length(wire_resistivities);
num_wire_flags = length(wire_material_flags);
num_scaling_factors = length(scaling_factors);
num_barrier_thicknesses = length(barrier_thicknesses);
num_barrier_resistivities = length(barrier_resistivities);
num_thermal_conductivities = length(thermal_conductivities);

total_configs = num_stacks * num_perms * num_thicks * num_freqs * num_cooling_configs * num_decaps * num_wire_resistivities * num_wire_flags * num_scaling_factors * num_barrier_thicknesses * num_barrier_resistivities * num_forced_powers * num_thermal_conductivities;

power = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);
power_density = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);
freq = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);
wire_power = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);
rep_power = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);
leakage_power = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);
dynamic_power = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);
temperature = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);
thickness = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);
npads = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);
cap_density = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);
Ltsv_m2 = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);
routable_design = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);

h_coeff = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);

num_metal_levels = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);
wire_pitch = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);

ild_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);
psn_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);
power_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);
wire_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);
chip_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);
tsv_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors,num_barrier_thicknesses,num_barrier_resistivities,num_forced_powers,num_thermal_conductivities);

%% Parameter sweep
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
                                    for bar_thick_ind = 1:num_barrier_thicknesses
                                        for bar_res_ind = 1:num_barrier_resistivities
                                            for forced_power_ind = 1:num_forced_powers
                                                for k_ind = 1:num_thermal_conductivities
                                                    cur_config = cur_config + 1;
                                                    die_thickness = thicknesses(thind);
                                                    num_layers_per_block = tiers(nind);
                                                    epsrd = rel_permittivities(pind);
                                                    fmax = frequencies(freq_ind);
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
                                                    fprintf('== bar thick: %d/%d \t=====\n',bar_thick_ind,num_barrier_thicknesses);
                                                    fprintf('==   bar res: %d/%d \t=====\n',bar_res_ind,num_barrier_resistivities);
                                                    fprintf('==       pow: %d/%d \t=====\n',forced_power_ind,num_forced_powers);
                                                    fprintf('==         K: %d/%d \t=====\n',k_ind,num_thermal_conductivities);
                                                    fprintf('==   Overall: %d/%d \t=====\n',cur_config,total_configs);
                                                    fprintf('===============================\n')

                                                    %% Rescale dimensions if necessary
                                                    compression_factor = scaling_factors(scaling_ind); % linear scaling factor. 1 = actual 32nm design, 4.57 = equivalent 7nm SB
                                                    Ach_mm2_scaled = Ach_mm2/compression_factor^2;
                                                    gate_pitch_scaled = gate_pitch/compression_factor;
                                                    min_pitch_scaled = min_pitch/compression_factor;
                                                    w_trans_scaled = w_trans/compression_factor;

                                                    % Let Vdd scale with scaling_ind
                                                    Vdd = xcm.get_nth_or_last( Vdd_vec, scaling_ind);

                                                    %% define parameters
                                                    [core.chip, core.transistor, core.gate, core.tsv, core.wire, core.psn, core.heat] = generate_basic_processor_settings(rent_exp,num_layers_per_block,Ng_core,Ach_mm2_scaled,gate_pitch_scaled,min_pitch_scaled,Vdd,fmax,w_trans_scaled);

                                                    %% Tweak wiring parameters
                                                    core.gate.output_resistance = 3e3/compression_factor; % Ohm
                                                    core.transistor.capacitance = 1e-15*3*w_trans_scaled*1e6; % ITRS projection is 1fF/um of gate width. This is an estimate for pMOS transistor capacitance

                                                    %core.wire.repeater_fraction = [0.3]; % 1 is default from gen_basic_proc_settings
                                                    %core.wire.routing_efficiency = [0.6]; % 0.4 is default from gen_basic_proc_settings
                                                    %core.wire.repeater_fraction = [0.4]; % 1 is default from gen_basic_proc_settings
                                                    %core.wire.routing_efficiency = [0.5]; % 0.4 is default from gen_basic_proc_settings
                                                    %core.wire.repeater_max_area_fraction = 0.2; % (-) Fraction of chip area that can be consumed by repeater/buffer gates
                                                    %core.wire.repeater_via_max_area_fraction = 0.05; % (-) Fraction of routable wire area that can be consumed by vias for repeater connections
                                                    %core.transistor.leakage_current_per_micron = 100e-9; %(A/um) % 32nm IOFF

                                                    wire_flag_key = bin2dec(wire_material_flags(wire_flag_ind));
                                                    alt_met_flag = bitget(wire_flag_key,1);
                                                    graphene_flag = bitget(wire_flag_key,2);
                                                    core.wire.use_graphene = graphene_flag;
                                                    core.wire.use_em_resistant_metal = alt_met_flag;   % (-) Allow or disallow use of electromigration-resistant metals below a specified minimum pitch
                                                    core.wire.min_non_em_width = 25e-9; % (m) If use_em_resistant_metal is set to 1, Cu resistivity will be replaced with wire.alt_resistivity_em below this pitch
                                                    core.wire.barrier_thickness = barrier_thicknesses(bar_thick_ind);
                                                    core.wire.barrier_resistivity = barrier_resistivities(bar_res_ind);
                                                    core.tsv.barrier_thickness = barrier_thicknesses(bar_thick_ind);
                                                    core.tsv.barrier_resistivity = barrier_resistivities(bar_res_ind);
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
                                                    core.psn.power_tsv_width = power_tsv_width;

                                                    if (simulation.force_power == 1)
                                                        core.chip.power_forced_val = power_forced_vec(forced_power_ind);
                                                    else
                                                        core.chip.power_forced_val = -1;
                                                    end


                                                    if (strcmp(cooling_configs{cind}, 'up') )
                                                        core.heat.up = heat_fluxes(cind);        % above chip
                                                        core.heat.down = core.heat.h_package;
                                                        core.heat.d = core.heat.h_package;
                                                    elseif (strcmp(cooling_configs{cind}, 'down') )
                                                        core.heat.up = core.heat.h_package;
                                                        core.heat.down = heat_fluxes(cind);
                                                        core.heat.d = core.heat.h_package;
                                                    elseif (strcmp(cooling_configs{cind}, 'down_all') )
                                                        core.heat.up = core.heat.h_package;
                                                        core.heat.down = heat_fluxes(cind);
                                                        core.heat.d = heat_fluxes(cind);
                                                    else % default to top-side heat sink
                                                        core.heat.up = heat_fluxes(cind);        % above chip
                                                        core.heat.down = core.heat.h_package;
                                                        core.heat.d = core.heat.h_package;
                                                    end

                                                    core.heat.k_underfill = thermal_conductivities(k_ind);
                                                    temperature_target = temperature_targets(cind);
                                                    core.chip.temperature = temperature_target;

                                                    %% If we're searching for frequencies below a certain temperature
                                                    if (simulation.freq_binsearch == 1)
                                                        core = find_thermally_limited_max_frequency(core, simulation, temperature_target);
                                                    elseif (simulation.power_binsearch == 1)
                                                        core_init = core;
                                                        simulation.ignore_leakage = 1;
                                                        core = find_thermally_limited_max_frequency(core, simulation, temperature_target);
                                                        power_target_W = core.power.total;
                                                        temperature_target_C = simulation.freq_binsearch_target;
                                                        fprintf('Beginning power search...\n')
                                                        simulation.ignore_leakage = 0;
                                                        core = design_for_power_target(power_target_W, temperature_target_C, core_init, simulation);
                                                    elseif (simulation.heat_transfer_binsearch == 1)
                                                        core = find_heat_transfer_coeff_for_target_temp(core, simulation);
                                                    else
                                                        %% calculate block parameters
                                                        [core.chip, core.power, core.tsv, core.wire, core.repeater, core.psn] = codesign_block(core.chip,core.tsv,core.gate,core.transistor,core.wire,core.heat,core.psn,simulation);
                                                    end
                                                    power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind) = core.power.total;
                                                    power_density(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind) = core.power.density;
                                                    freq(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind) = 1/core.chip.clock_period;

                                                    wire_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind) = core.power.wiring;
                                                    rep_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind) = core.power.repeater;
                                                    dynamic_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind) = core.power.repeater_dynamic + core.power.dynamic;
                                                    leakage_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind) = core.power.repeater_leakage + core.power.leakage;
                                                    temperature(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind) = core.chip.temperature;
                                                    thickness(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind) = core.chip.thickness;
                                                    npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind) = core.psn.Npads;
                                                    h_coeff(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind) = core.heat.up;

                                                    num_metal_levels(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind) = length(core.wire.pn);
                                                    wire_pitch{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind} = core.wire.pn;

                                                    ild_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind} = core.chip.iidf;
                                                    psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind} = core.psn;
                                                    power_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind} = core.power;
                                                    wire_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind} = core.wire;
                                                    chip_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind} = core.chip;
                                                    tsv_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind} = core.tsv;

                                                    routable_design(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind) = core.wire.routable;

                                                    % Get L and C density for the power supply network
                                                    if ( (simulation.skip_psn_loops == 0) && (core.wire.routable == 1) )
                                                        psn_inductance_density_m2 = psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind}.Ltsv/psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind}.l_unit_cell^2;
                                                        psn_capacitance_density_m2 = psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind}.cap_density;
                                                    else
                                                        psn_inductance_density_m2 = 0;
                                                        psn_capacitance_density_m2 = 0;
                                                    end

                                                    Ltsv_m2(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind) = psn_inductance_density_m2;
                                                    cap_density(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind) = psn_capacitance_density_m2;
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
        end
    end
end


t_sweep_stop = cputime;
fprintf('\nTotal time elapsed for parameter sweep: %.3g seconds\t(%.3g minutes)\n\n',(t_sweep_stop-t_sweep_start),(t_sweep_stop-t_sweep_start)/60);
