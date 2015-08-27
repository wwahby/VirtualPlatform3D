function core = find_thermally_limited_max_frequency(core, simulation)

max_gens = simulation.freq_binsearch_max_gens;
target_max_value = simulation.freq_binsearch_target;
target_cur_value = target_max_value*1; % start with something invalid so we run it at least once
abs_err = abs(target_max_value - target_cur_value);
tolerance = simulation.freq_binsearch_raw_tol;

%% Find bounds
% Start with initial guess g0
% if T(g(i)) < target, g(i+1) = A*g(i)
% if T(g(i)) > target, g(i+1) = g(i)/A
% Stop when T(g(i)) flips relation with target (i.e. is less than when it
% was previously greater than, etc

A = 10; % amount to multiply or divide by at each generation

keep_going = true;
within_tol = false;
cur_freq = simulation.freq_binsearch_initial_guess;

gen_init = 0;
dir_sign = false;
prev_dir_sign = false;
prev_freq = cur_freq;
prefactor = 1;

time_init_start = cputime;
while(keep_going)
    core.chip.clock_period = 1/cur_freq;
    [core.chip, core.power, core.tsv, core.wire, core.repeater, core.psn] = codesign_block(core.chip,core.tsv,core.gate,core.transistor,core.wire,core.heat,core.psn,simulation);
    cur_temp = core.chip.temperature;
    abs_err = abs(target_max_value - cur_temp);
    within_tol = (abs_err <= tolerance);
    
    prev_dir_sign = dir_sign;
    if (cur_temp < target_max_value)
        prefactor = A;
        dir_sign = true;
    elseif (cur_temp > target_max_value)
        prefactor = 1/A;
        dir_sign = false;
    else % cur_temp == target_max_value -- shouldn't happen, but you never know
        prefactor = 1;
        dir_sign = false;
        % [ FIX ] need to break here so we don't get stuck
    end
    
    if ( gen_init > 0) % the second run is the first one where we'll know if we can stop
        if ( dir_sign ~= prev_dir_sign)
            keep_going = false;
            min_bound = min(prev_freq, cur_freq);
            max_bound = max(prev_freq, cur_freq);
        end
    end
    
    keep_going = keep_going && (~within_tol); % also stop if we're within the tolerance value;
    fprintf('Initial Search Gen %d: \t Freq: %.3g \t Temp: %.4d\n\n', gen_init, cur_freq, cur_temp);
    if (keep_going)
        prev_freq = cur_freq;
        cur_freq = prefactor*prev_freq;
        gen_init = gen_init + 1;
    end
end

% At this point min_bound and max_bound should be set
% We may have accidentally hit on something close enough to accept, in that
% case just stop
% Otherwise binary search with the bounds we found.
if (~within_tol)
    %% Binary search
    left = min_bound;
    right = max_bound;
    mid = 1/2*(left+right);
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
        fprintf('Binsrch Gen %d: \t Freq: %.3g \t Temp: %.4d\n\n',num_gens, mid, target_cur_value);
    end
end
time_bin_stop = cputime;