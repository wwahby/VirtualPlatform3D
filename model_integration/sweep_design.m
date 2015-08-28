function sweep_data = sweep_design(design, sweep, simulation)
t_sweep_start = cputime;
%% Unpack design settings
rent_exp = design.rent_exp;
Ng_core = design.Ng_core;
Ach_mm2 = design.Ach_mm2;
gate_pitch = design.gate_pitch;
min_pitch = design.min_pitch;
Vdd = design.Vdd;
w_trans = design.w_trans;

%% Unpack sweep settings
tiers = sweep.tiers;
thicknesses = sweep.thicknesses;
force_thickness = sweep.force_thickness;
rel_permittivities = sweep.rel_permittivities;
frequencies = sweep.frequencies;
heat_fluxes = sweep.heat_fluxes;
decap_ratios = sweep.decap_ratios;
wire_resistivities = sweep.wire_resistivities;
wire_material_flags = sweep.wire_material_flags;
scaling_factors = sweep.scaling_factors;

%% Initialize sweep variables
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
num_scaling_factors = length(scaling_factors);
total_configs = num_stacks * num_perms * num_thicks * num_freqs * num_cooling_configs * num_decaps * num_wire_resistivities * num_wire_flags * num_scaling_factors;

sweep_data.power = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
sweep_data.power_density = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
sweep_data.freq = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
sweep_data.wire_power = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
sweep_data.rep_power = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
sweep_data.temp = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
sweep_data.thickness = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
sweep_data.npads = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
sweep_data.cap_density = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
sweep_data.Ltsv_m2 = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);

sweep_data.num_metal_levels = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
sweep_data.wire_pitch = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);

sweep_data.ild_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
sweep_data.psn_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
sweep_data.power_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
sweep_data.wire_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
sweep_data.chip_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);
sweep_data.tsv_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs,num_wire_resistivities,num_wire_flags,num_scaling_factors);

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
                                    fprintf('==   Overall: %d/%d \t=====\n',cur_config,total_configs);
                                    fprintf('===============================\n')
                                    
                                    %% Rescale dimensions if necessary
                                    compression_factor = scaling_factors(scaling_ind); % linear scaling factor. 1 = actual 32nm design, 4.57 = equivalent 7nm SB
                                    Ach_mm2_scaled = Ach_mm2/compression_factor^2;
                                    gate_pitch_scaled = gate_pitch/compression_factor;
                                    min_pitch_scaled = min_pitch/compression_factor;
                                    w_trans_scaled = w_trans/compression_factor;
                                    
                                    %% define parameters
                                    [core.chip, core.transistor, core.gate, core.tsv, core.wire, core.psn, core.heat] = generate_basic_processor_settings(rent_exp,num_layers_per_block,Ng_core,Ach_mm2_scaled,gate_pitch_scaled,min_pitch_scaled,Vdd,fmax,w_trans_scaled);

                                    %% Tweak wiring parameters
                                    core.gate.output_resistance = 8e3*compression_factor; % Ohm
                                    core.transistor.capacitance = 1e-15*1e6*3*w_trans; % ITRS projection is 1fF/um of gate width. This is an estimate for pMOS transistor capacitance

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
                                        core = find_thermally_limited_max_frequency(core, simulation);
                                    else
                                        %% calculate block parameters
                                        [core.chip, core.power, core.tsv, core.wire, core.repeater, core.psn] = codesign_block(core.chip,core.tsv,core.gate,core.transistor,core.wire,core.heat,core.psn,simulation);
                                    end
                                    sweep_data.power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.power.total;
                                    sweep_data.power_density(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.power.density;
                                    sweep_data.freq(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = 1/core.chip.clock_period;
                                    

                                    sweep_data.wire_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.power.wiring;
                                    sweep_data.rep_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.power.repeater;
                                    sweep_data.temp(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.chip.temperature;
                                    sweep_data.thickness(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.chip.thickness;
                                    sweep_data.npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.psn.Npads;
                                    
                                    sweep_data.num_metal_levels(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = length(core.wire.pn);
                                    sweep_data.wire_pitch{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind} = core.wire.pn;

                                    sweep_data.ild_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind} = core.chip.iidf;
                                    sweep_data.psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind} = core.psn;
                                    sweep_data.power_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind} = core.power;
                                    sweep_data.wire_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind} = core.wire;
                                    sweep_data.chip_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind} = core.chip;
                                    sweep_data.tsv_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind} = core.tsv;

                                    if (simulation.skip_psn_loops == 0)
                                        sweep_data.Ltsv_m2(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind}.Ltsv/psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind}.l_unit_cell^2;
                                        sweep_data.cap_density(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind}.cap_density;
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
