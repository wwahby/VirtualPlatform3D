function [chip, power, tsv, wire, repeater, psn] = codesign_block(chip,tsv,gate,transistor,wire,heat,psn,simulation)
time_start = cputime;
% Power and signal codesign
%% constants
eps0 = 8.854e-12; % (F/m) vacuum permittivity
mu0 = 4*pi*1e-7; % (H/m) Vacuum permeability

%% Initializations
temperature_previous = chip.temperature;
temperature_current = chip.temperature;
temperature_converged = 0;
temperature_iterations = 1;
temperature_iterations_max = 10;
temperature_change_tolerance = 1;

%% TSV number determination
% Inputs:
%   Rent parameters
%   Number of logic gates
%   Number of layers

% disp(' ')
% disp('Estimating TSV requirements...')
[nt_max, nt_tot, nt_to, nt_through, Tacmat] = xcm.estimate_tsvs_required(chip.num_gates,chip.num_layers,chip.rent_k,chip.rent_p,chip.alpha);

%% TSV Sizing
% Inputs:
%   Area available
%   Max area for TSVs
%   TSV aspect ratio
% disp('Sizing TSVs...')
[w_tsv_m, h_tsv_m] = xcm.size_tsvs(chip.area_total/chip.num_layers, tsv.max_area_fraction, nt_max, tsv.aspect_ratio );
%h_tsv_gp = round(h_tsv_gp);

if (simulation.force_thickness == 1)
    h_tsv_m = chip.thickness_nominal;
    w_tsv_m = h_tsv_m/tsv.aspect_ratio;
end
h_tsv_gp = ceil(h_tsv_m/chip.gate_pitch);
w_tsv_gp = ceil(w_tsv_m/chip.gate_pitch);

tsv.width_m = w_tsv_m;
tsv.height_m = h_tsv_m;
tsv.width_gp = w_tsv_gp;
tsv.height_gp = h_tsv_gp;
tsv.per_layer = nt_tot;
tsv.max_tsvs_per_layer = nt_max;

if(tsv.height_m > 0)
    chip.thickness = tsv.height_m;
else
    chip.thickness = chip.thickness_nominal;
end

%% Check if we need to run/rerun the interconnect and thermal modules
initial_run = (temperature_iterations == 1);
temp_below_sanity_limit = (chip.temperature < simulation.insanity_temperature);
rerun_temperature = (~temperature_converged) && (temperature_iterations <= temperature_iterations_max) && (simulation.iterate_temperature) && (~simulation.skip_thermal) && (temp_below_sanity_limit);

run_xcm_and_thermal_modules = (initial_run || rerun_temperature);

while(run_xcm_and_thermal_modules)
    
    temperature_previous = temperature_current;
    %% System determination
    % Run WLD + WLA + RI to get power estimate
    % disp('Generating system...')
    [chip, power, wire, repeater, tsv] = xcm.gen_design(chip,tsv,gate,transistor,wire,simulation);

    %Pdens = power.total/chip.area_per_layer_m2;
    %power.density = power.total/chip.area_total;

    %% Thermal module -- Find actual system temperature

    if (simulation.force_power == 1)
        % Find frequency which will give us correct power
        Ptot_new = chip.power_forced_val;

        Plk_tot = power.leakage + power.repeater_leakage;
        Pdyn_tot = power.dynamic + power.repeater_dynamic + power.wiring;

        cur_freq = 1/chip.clock_period;
        new_freq = (Ptot_new - Plk_tot)/(Pdyn_tot/cur_freq);
        chip.clock_period = 1/new_freq;

        Pdyn_logic = power.dynamic;
        power.dynamic = Pdyn_logic/cur_freq*new_freq;
        power.repeater_dynamic = power.repeater_dynamic/cur_freq*new_freq;
        power.wiring = power.wiring/cur_freq*new_freq;

        power.repeater = power.repeater_leakage + power.repeater_dynamic;
        power.total = power.repeater + power.dynamic + power.leakage + power.wiring;
        power.density = power.total/chip.area_per_layer_m2;
    end
    power_per_layer = power.total/chip.num_layers;
    power_therm_vec = ones(1,chip.num_layers)*power_per_layer;  %power dissipation of each die
    layer_area = chip.area_per_layer_m2;
    side_length = sqrt(layer_area);
    chip_width = side_length;
    chip_height = side_length;

    package_width = 3.5*chip_width;
    package_height = 3.5*chip_height;

    %power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....
    map_row = [0     0     side_length    side_length     power_per_layer];
    map = zeros(chip.num_layers,5);
    blk_num = zeros(chip.num_layers,1);
    for i =1:chip.num_layers
        map(i,:) = map_row;
        blk_num(i,1) = 1;
    end
    %blk_num is for splitting the power maps of each die

    if(~simulation.skip_thermal)
        [max_temp temp_vec] = get_stack_temperature(chip.num_layers,chip.thickness,wire,tsv,chip_width,chip_height,package_width,package_height,heat,simulation,map,blk_num,power_therm_vec);

        chip.temperature_vec = temp_vec;
        chip.temperature = max_temp;

        temperature_current = chip.temperature;
        temperature_change = temperature_current - temperature_previous;
        fprintf('\tPrevious Temperature: %.5g\tCurrent Temperature: %.5g\n',temperature_previous, temperature_current);

        temperature_converged = ( abs(temperature_change) < temperature_change_tolerance);
    else
        fprintf('\tThermal analysis SKIPPED due to simulation flag!\n')
    end
    temperature_iterations = temperature_iterations + 1;
    
    initial_run = (temperature_iterations == 1);
    temp_below_sanity_limit = (chip.temperature < simulation.insanity_temperature);
    rerun_temperature = (~temperature_converged) && (temperature_iterations <= temperature_iterations_max) && (simulation.iterate_temperature) && (~simulation.skip_thermal) && (temp_below_sanity_limit) && (simulation.iterate_temperature);
    run_xcm_and_thermal_modules = (initial_run || rerun_temperature);
    chip.temperature_exceeds_sanity_limit = ~temp_below_sanity_limit;
end



%% Power Supply Network Determination
if (wire.routable == 0)
    fprintf('\tPSN analysis SKIPPED due to unroutable design!\n')
elseif (~simulation.skip_psn_loops)
    psn = determine_power_tsv_requirements(tsv,psn,power,wire,chip);
else
    fprintf('\tPSN analysis SKIPPED due to simulation flag!\n')
end

temperature_K = chip.temperature + 273.15;
decap_frac_per_die = psn.decap_area_fraction * ones(1,chip.num_layers);
power_per_layer = power.total/chip.num_layers;
power_per_die = ones(1,chip.num_layers)*power_per_layer;  %power dissipation of each die
resistivity_bulk = wire.resistivity;
if (simulation.run_transient_psn)
    [max_noise, max_noise_time, time_mat, voltage_mat] = pdnt.run_transient_power_analysis(chip.Vdd, chip.thickness, chip_width, chip_height, tsv, resistivity_bulk, power_per_die, decap_frac_per_die, temperature_K, heat );
    psn.transient.max_noise = max_noise;
    psn.transient.max_noise_time = max_noise_time;
    psn.transient.time_mat = time_mat;
    psn.transient.voltage_mat = voltage_mat;
end



%% Final report

time_stop = cputime;
time_elapsed = time_stop - time_start;
% disp(' ')
% disp('Final block parameters:')
% repstr = sprintf('\tNg_nom %d \t Ng_act: %d \t Atsv_nom: %.3g \t Atsv_act: %.3g \n\tN_tsvs: %d \t Npads_pow %d \t psn_nom %.4g \t psn_act %.4g', ...
%                   chip.num_gates, chip.Ng_actual, tsv.max_area_fraction, tsv.actual_area_fraction, tsv.num, psn.Npads, psn.noise_target, psn_max);
% disp(repstr)
% repstr = sprintf('\th_tsv_um: %.4g \t w_tsv_um: %.4g',h_tsv_m/1e-6,w_tsv_m/1e-6);
% disp(repstr);

fprintf('\tTotal time elapsed for block design: %.3g seconds\n',time_elapsed)

