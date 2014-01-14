function [chip power tsv wire repeater psn] = codesign_system(chip,tsv,gate,transistor,wire,psn,simulation)
% Power and signal codesign
%% constants
eps0 = 8.854e-12; % (F/m) vacuum permittivity
mu0 = 4*pi*1e-7; % (H/m) Vacuum permeability

%% TSV number determination
% Inputs:
%   Rent parameters
%   Number of logic gates
%   Number of layers
disp(' ')
disp('Estimating TSV requirements...')
[nt_max nt_tot nt_to nt_through Tacmat] = estimate_tsvs_required(chip.num_gates,chip.num_layers,chip.rent_k,chip.rent_p,chip.alpha);

%% TSV Sizing
% Inputs:
%   Area available
%   Max area for TSVs
%   TSV aspect ratio
disp('Sizing TSVs...')
[w_tsv_gp h_tsv_gp] = size_tsvs(chip.area_total/chip.num_layers, tsv.max_area_fraction, nt_max, tsv.aspect_ratio );
%h_tsv_gp = round(h_tsv_gp);
h_tsv_m = h_tsv_gp*chip.gate_pitch;
w_tsv_m = w_tsv_gp*chip.gate_pitch;

tsv.width_m = w_tsv_m;
tsv.height_m = h_tsv_m;
tsv.width_gp = w_tsv_gp;
tsv.height_gp = h_tsv_gp;

%% System determination
% Run WLD + WLA + RI to get power estimate
disp('Generating system...')
%[ iidf l Ln pn pn_orig Cxc Ltot Cn Pdyn Plk Pw Prep Ng_act N_tsvs T_tsvs Atf_act ] = gen_design_old(Ng,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);
[chip power wire repeater tsv] = gen_design(chip,tsv,gate,transistor,wire,simulation);

%Pdens = power.total/chip.area_per_layer_m2;
power.density = power.total/chip.area_total;
psn.pitch_tsv = tsv.pitch_m*10; % [FIX] Power TSVs aren't going to be on the same pitch as signal TSVs

%% System temperature check may go here

%% Power noise estimation
% Inputs:
%   Total power (from previous step)
%   TSV geometry from TSV Sizing step
disp('Evaluating power supply noise...')

psn_iterations = 1;
rho_m = wire.resistivity(1); % use top layer
mu_m = wire.permeability_rel*mu0;
psn_max = calc_psn(psn,power,chip,tsv,rho_m,mu_m,chip.temperature);

dispstr = sprintf('\tpsn_runs: %d\tNpads: %d\tpsn_target: %d\tpsn_max: %d',psn_iterations,psn.Npads, psn.noise_target, psn_max);
disp(dispstr)


%% Home in on the best number of power TSVs to use to meet the target
psn.noise = psn_max;
mismatch_norm = psn.noise/psn.noise_target;
psn_iterations_max = 20;
psn_mismatch_tolerance = 0.05; % 5% mismatch is OK
while ( (abs(mismatch_norm-1) > psn_mismatch_tolerance) && (psn_iterations < psn_iterations_max) )
    % more pads -> less noise, and vice versa
    psn.Npads_1d = round(sqrt(psn.Npads*mismatch_norm));
    psn.Npads = psn.Npads_1d^2;
    psn_max = calc_psn(psn,power,chip,tsv,rho_m,mu_m,chip.temperature);
    psn.noise = psn_max;
    mismatch_norm = psn.noise/psn.noise_target;
    dispstr = sprintf('\tpsn_runs: %d\tNpads: %d\tpsn_target: %d\tpsn_max: %d\tmismatch_norm: %.3g',psn_iterations,psn.Npads, psn.noise_target, psn_max,mismatch_norm);
	disp(dispstr)
end

%% Final report
disp(' ')
disp('Final system parameters:')
repstr = sprintf('\tNg_nom %d \t Ng_act: %d \t Atsv_nom: %.3g \t Atsv_act: %.3g \n\tN_tsvs: %d \t Npads_pow %d \t psn_nom %.4g \t psn_act %.4g', ...
                  chip.num_gates, chip.Ng_actual, tsv.max_area_fraction, tsv.actual_area_fraction, tsv.num, psn.Npads, psn.noise_target, psn_max);
disp(repstr)
repstr = sprintf('\th_tsv_um: %.4g \t w_tsv_um: %.4g',h_tsv_m/1e-6,w_tsv_m/1e-6);
disp(repstr);

