function psn = determine_power_tsv_requirements(tsv,psn,power,wire,chip)
%% Power noise estimation
% Inputs:
%   Total power (from previous step)
%   TSV geometry from TSV Sizing step

% constants
eps0 = 8.854e-12; % (F/m) vacuum permittivity
mu0 = 4*pi*1e-7; % (H/m) Vacuum permeability

% disp('Evaluating power supply noise...')
%psn.pitch_tsv = tsv.pitch_m*10; % [FIX] Power TSVs aren't going to be on the same pitch as signal TSVs
psn_mismatch_tolerance = psn.mismatch_tolerance; % 5% mismatch is OK
psn_iterations = 1;

rho_m = wire.rho_vec(end); % use top layer
mu_m = wire.permeability_rel*mu0;

psn.noise_target = psn.noise_fraction * chip.Vdd;

%% Run twice to jump close to where we need to be
norm_err_min = 1e3; % setting this to an arbitrarily high value to start

% fprintf('\tLimiting search space...\n')
[psn_max, RTSV, LTSV, cap_density, l_unit_cell] = power_noise.calc_psn(psn,power,chip,tsv,rho_m,mu_m,chip.temperature);
psn.noise = psn_max;
psn.Rtsv = RTSV;
psn.Ltsv = LTSV;
psn.cap_density = cap_density;
psn.l_unit_cell = l_unit_cell;
psn.noise = psn_max;
mismatch_norm = psn.noise/psn.noise_target;
psn_target = psn.noise_target;
rel_err = (psn_max - psn_target)/psn_target;
norm_err = abs(rel_err);

if (norm_err < norm_err_min) % new high score!
    norm_err_min = norm_err;
    psn_best = psn;
end

% dispstr = sprintf('\tpsn_runs: %d\tNpads: %d\tpsn_target: %d\tpsn_max: %d\trel_err: %.3g',0,psn.Npads, psn.noise_target, psn_max,rel_err);	
% disp(dispstr)


new_pads_1d = round(sqrt(psn.Npads*mismatch_norm));
psn.Npads_1d = new_pads_1d;
psn.Npads = psn.Npads_1d^2;
[psn_max, RTSV, LTSV, cap_density, l_unit_cell] = power_noise.calc_psn(psn,power,chip,tsv,rho_m,mu_m,chip.temperature);
psn.noise = psn_max;
psn.Rtsv = RTSV;
psn.Ltsv = LTSV;
psn.cap_density = cap_density;
psn.l_unit_cell = l_unit_cell;
rel_err = (psn_max - psn_target)/psn_target;
norm_err = abs(rel_err);

if (norm_err < norm_err_min) % new high score!
    norm_err_min = norm_err;
    psn_best = psn;
end

npads_first_linear_bound = psn.Npads_1d;
err_first_linear_bound = norm_err;
psn_first_linear_bound = psn;

% dispstr = sprintf('\tpsn_runs: %d\tNpads: %d\tpsn_target: %.3g\tpsn_max: %.3g\trel_err: %.3g',0,psn.Npads, psn.noise_target, psn_max,rel_err);	
% disp(dispstr)


%% Find bounds for binary search
% fprintf('\tFinding bounds for binary search...\n')
psn_target = psn.noise_target;
npads_1d = psn.Npads_1d;
npads = psn.Npads;
if (psn_max < psn_target)
    while((psn_max < psn_target) && (npads_1d > 1));
        npads_1d_old = npads_1d;
        npads_old = npads;
        npads_1d = round(1/2*npads_1d);
        psn.Npads_1d = npads_1d;
        psn.Npads = psn.Npads_1d^2;
        npads = psn.Npads;
        [psn_max RTSV LTSV cap_density l_unit_cell] = power_noise.calc_psn(psn,power,chip,tsv,rho_m,mu_m,chip.temperature);
        psn.noise = psn_max;
        psn.Rtsv = RTSV;
        psn.Ltsv = LTSV;
        psn.cap_density = cap_density;
        psn.l_unit_cell = l_unit_cell;

        rel_err = (psn_max - psn_target)/psn_target;
        norm_err = abs(rel_err);

        if (norm_err < norm_err_min) % new high score!
            norm_err_min = norm_err;
            psn_best = psn;
        end

%         fprintf('\tpsn_runs: %d\tNpads: %d\tpsn_target: %.3g\tpsn_max: %.3g\trel_err: %.3g\n',0,psn.Npads, psn.noise_target, psn_max,rel_err);
    end

%     lbnd = npads_1d;
%     rbnd = npads_1d_old;
    lbnd = npads;
    rbnd = npads_old;
elseif(psn_max > psn_target)
    while(psn_max > psn_target)

        npads_1d_old = npads_1d;
        npads_old = npads;
        npads_1d = 2*npads_1d;
        psn.Npads_1d = npads_1d;
        psn.Npads = psn.Npads_1d^2;
        npads = psn.Npads;
        [psn_max RTSV LTSV cap_density l_unit_cell] = power_noise.calc_psn(psn,power,chip,tsv,rho_m,mu_m,chip.temperature);
        psn.noise = psn_max;
        psn.Rtsv = RTSV;
        psn.Ltsv = LTSV;
        psn.cap_density = cap_density;
        psn.l_unit_cell = l_unit_cell;

        rel_err = (psn_max - psn_target)/psn_target;
        norm_err = abs(rel_err);

        if (norm_err < norm_err_min) % new high score!
            norm_err_min = norm_err;
            psn_best = psn;
        end

%         fprintf('\tpsn_runs: %d\tNpads: %d\tpsn_target: %.3g\tpsn_max: %.3g\trel_err: %.3g\n',0,psn.Npads, psn.noise_target, psn_max,rel_err);
    end
%     lbnd = npads_1d_old;
%     rbnd = npads_1d;
    lbnd = npads_old;
    rbnd = npads;
end

% npads_second_linear_bound = psn.Npads_1d;
% err_second_linear_bound = norm_err;
% psn_second_linear_bound = psn;

%% Binary search

gen_ind = 1;
max_gens = 16;
tol = psn_mismatch_tolerance;
rel_err = (psn_max - psn_target)/psn_target;
norm_err = abs(rel_err);
while ((norm_err > tol) && (gen_ind < max_gens) && (rbnd > lbnd+1))

    mid = round((lbnd+rbnd)/2);
    %psn.Npads_1d = mid;
    psn.Npads = mid;
    %psn.Npads = psn.Npads_1d^2;
    [psn_max, RTSV, LTSV, cap_density, l_unit_cell, output_cell] = power_noise.calc_psn(psn,power,chip,tsv,rho_m,mu_m,chip.temperature);
    psn.noise = psn_max;
    psn.Rtsv = RTSV;
    psn.Ltsv = LTSV;
    psn.cap_density = cap_density;
    psn.l_unit_cell = l_unit_cell;
    psn.output_cell = output_cell;

    if (psn_max > psn_target) % need more tsvs
        lbnd = mid;
    else % psn_max < psn_target, can get away with less tsvs
        rbnd = mid;
    end

    rel_err = (psn_max - psn_target)/psn_target ;
    norm_err = abs(rel_err);

    if (norm_err < norm_err_min) % new high score!
        norm_err_min = norm_err;
        psn_best = psn;
    end

%     fprintf('\tpsn_runs: %d\tNpads: %d\tpsn_target: %d\tpsn_max: %d\trel_err: %.3g\n',gen_ind,psn.Npads, psn.noise_target, psn_max,rel_err);
    gen_ind = gen_ind + 1;
end

if (psn.power_tsv_width == -1)
    psn.power_tsv_width = tsv.width_m;
end


 fprintf('\tPSN Done! \tNpads: %d\tpsn_target: %d\tpsn_max: %d\tnorm_err: %.3g\n',psn.Npads, psn.noise_target, psn_max,norm_err_min);
