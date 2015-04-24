function [chip, power, tsv, wire, repeater, psn] = codesign_block(chip,tsv,gate,transistor,wire,heat,psn,simulation)
time_start = cputime;
% Power and signal codesign
%% constants
eps0 = 8.854e-12; % (F/m) vacuum permittivity
mu0 = 4*pi*1e-7; % (H/m) Vacuum permeability

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
%% System determination
% Run WLD + WLA + RI to get power estimate
% disp('Generating system...')
[chip, power, wire, repeater, tsv] = xcm.gen_design(chip,tsv,gate,transistor,wire,simulation);

%Pdens = power.total/chip.area_per_layer_m2;
%power.density = power.total/chip.area_total;

%% Thermal module -- Find actual system temperature

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


[max_temp temp_vec] = get_stack_temperature(chip.num_layers,chip.thickness,wire,tsv,chip_width,chip_height,package_width,package_height,heat,simulation,map,blk_num,power_therm_vec);

chip.temperature_vec = temp_vec;
chip.temperature = max_temp;



%% Power Supply Network Determination
if (~simulation.skip_psn_loops)
    psn = determine_power_tsv_requirements(tsv,psn,power,wire,chip);
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

