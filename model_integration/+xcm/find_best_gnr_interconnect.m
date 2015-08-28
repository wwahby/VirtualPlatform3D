function [use_gnr, gnr_width, gnr_pitch, gnr_delay, R, C, C_pul] = ...
    find_best_gnr_interconnect( num_layers, gnr_length, delay_max, min_pitch, ...
                                width_fraction, pitch_orig, temp_K, mfp_defect, ...
                                rho_interlayer, prob_backscattering, Ef, contact_resistance, ...
                                epsrd, height_dielectric)
                            
% initialize binary search parameters
max_gens = 32; % maximum number of binary search generations to run
err = inf;
err_rel = inf;
rel_err_tolerance = 1e-2;
A = 10; % prefactor for finding intial binary search bounds

gnr_width = pitch_orig*width_fraction;
gnr_space = pitch_orig - gnr_width;
min_width = width_fraction*min_pitch;

last_valid_width = -1; % start with invalid width


%% First find bounds for binary search
% Start with guess and either go up or down by factor of A (10X typ) until
% we have bounds

keep_going = true;
within_tol = false;
cur_width = gnr_width;
cur_pitch = cur_width/width_fraction;
cur_space = cur_pitch - cur_width;

gen_init = 0;
dir_sign = false; % false = decrease, true = increase
prev_dir_sign = false;
prev_width = cur_width;
prefactor = 1;
cur_delay = -1;

while (keep_going)
    [delay_top, delay_side, R_top, R_top_alt, R_side, L_vec, C_gnr_vec, C_gnr_raw_vec, Nch_vec, mfp_eff_vec ] = ...
        xcm.calc_gnr_params_combined_multiple_widths( ...
                num_layers, cur_width, cur_space, gnr_length, ...
                temp_K, mfp_defect, rho_interlayer, prob_backscattering, ...
                Ef, contact_resistance, epsrd, height_dielectric );

    cur_delay = delay_top;
    
    err_raw = cur_delay - delay_max;
    err_rel = err_raw/delay_max;
    
    within_tol = ( abs(err_rel) <= rel_err_tolerance );
    fprintf('keep_going: %d \t within_tol: %d\n', keep_going, within_tol)

    prev_dir_sign = dir_sign;
    if (cur_delay < delay_max) % XC too fast, decrease size of wire
        prefactor = 1/A;
        dir_sign = false;
    elseif (cur_delay > delay_max) % XC too slow, increase size of wire
        prefactor = A;
        dir_sign = true;
    else % cur_temp == target_max_value -- shouldn't happen, but you never know
        prefactor = 1;
        dir_sign = false;
        keep_going = false; % Break out of this loop
    end
    
    if ( gen_init > 0) % the second run is the first one where we'll know if we can stop
        if ( dir_sign ~= prev_dir_sign)
            keep_going = false;
            min_bound = min(prev_width, cur_width);
            max_bound = max(prev_width, cur_width);
        end
    end
    
    keep_going = keep_going && (~within_tol); % also stop if we're within the tolerance value;
    fprintf('   Initial Search Gen %d: \t Width: %.3g \t Delay: %.3g \t Delay target: %.3g\n', gen_init, cur_width, cur_delay, delay_max);
    if (keep_going)
        prev_width = cur_width;
        
        cur_width = prefactor*prev_width;
        cur_pitch = cur_width/width_fraction;
        cur_space = cur_pitch - cur_width;
        
        gen_init = gen_init + 1;
    end
end

fprintf('out of the loop\n')

% At this point min_bound and max_bound should be set
% We may have accidentally hit on something close enough to accept, in that
% case just stop
% Otherwise binary search with the bounds we found.
if (~within_tol)
    %% Binary search
    left_width = min_bound;
    right_width = max_bound;
    mid_width = 1/2*(left_width+right_width);
    mid_pitch = mid_width/width_fraction;
    mid_space = mid_pitch - mid_width;
    
    num_gens = 0;
    time_bin_start = cputime;
    while ( (abs(err_rel) > rel_err_tolerance) && (num_gens < max_gens) )
        
        mid_width = 1/2*(left_width+right_width);
        mid_pitch = mid_width/width_fraction;
        mid_space = mid_pitch - mid_width;
        
        [delay_top, delay_side, R_top, R_top_alt, R_side, L_vec, C_gnr_vec, C_gnr_raw_vec, Nch_vec, mfp_eff_vec ] = ...
        xcm.calc_gnr_params_combined_multiple_widths( ...
                num_layers, mid_width, mid_space, gnr_length, ...
                temp_K, mfp_defect, rho_interlayer, prob_backscattering, ...
                Ef, contact_resistance, epsrd, height_dielectric );

        cur_delay = delay_top;
        err_raw = cur_delay - delay_max;
        err_rel = err_raw/delay_max;

        if (cur_delay > delay_max) % XC too slow, need to increase size
            left_width = mid_width;
        elseif (cur_delay < target_max_value) % XC too fast, need to decrease size
            right_width = mid_width ;
        else % XC exactly right, just decrease search interval. This isn't important as we'll just exit out of this loop in a moment
            left_width = 1/2*(left_width + mid_width);
            right_width = 1/2*(right_width + mid_width);
        end

        num_gens = num_gens + 1;
        fprintf('   Binsrch Gen %d: \t Freq: %.3g \t Temp: %.4d\n',num_gens, mid, target_cur_value);
    end
end

time_bin_stop = cputime;
time_binsearch = time_bin_stop - time_bin_start;
fprintf('   GNR Insertion Done! Time elapsed: %.3g', time_binsearch);

gnr_width = mid_width;
gnr_pitch = gnr_width/width_fraction;
gnr_delay = cur_delay;
R = R_top;
C = C_gnr_vec;
C_pul = C_gnr_vec/gnr_length;




%% Old
% 
% % First check that GNRs can meet delay constraint at original size
% [delay_top, delay_side, R_top, R_top_alt, R_side, L_vec, C_gnr_vec, C_gnr_raw_vec, Nch_vec, mfp_eff_vec ] = ...
%     xcm.calc_gnr_params_combined_multiple_widths( ...
%         num_layers, gnr_width, gnr_space, gnr_length, temp_K, mfp_defect, rho_interlayer, prob_backscattering, ...
%         Ef, contact_resistance, epsrd, height_dielectric );
% 
% 
% if (delay_top > delay_max) % GNR no better than Cu -- don't bother trying to shrink the wires
%     use_gnr = 0;
%     best_R = R_top;
%     best_C = C_gnr_vec;
%     best_C_pul = C_gnr_vec/gnr_length;
%     best_delay = delay_top;
% else
%     use_gnr = 1;
%     
%     % if min pitch works, just end here and move along
%     [delay_top delay_side R_top R_top_alt R_side L_vec C_gnr_vec C_gnr_raw_vec Nch_vec mfp_eff_vec ] = ...
%     xcm.calc_gnr_params_combined_multiple_widths( ...
%         num_layers, min_width, gnr_space, gnr_length, temp_K, mfp_defect, rho_interlayer, prob_backscattering, ...
%         Ef, contact_resistance, epsrd, height_dielectric );
%     delay_gnr = delay_top;
%     pass = (delay_gnr <= delay_max);
%     
%     if(pass)
%         last_valid_width = min_width;
% 
%         best_R = R_top;
%         best_C = C_gnr_vec;
%         best_C_pul = C_gnr_vec/gnr_length;
%         best_delay = delay_gnr;
%     else % check the rest of the widths
% 
%         last_valid_width = gnr_width;
%         num_gens = 1; % binary search generation counter
% 
%         best_R = R_top;
%         best_C = C_gnr_vec;
%         best_C_pul = C_gnr_vec/gnr_length;
%         best_delay = delay_top;
% 
%         lbnd = min_width;
%         rbnd = gnr_width;
%         mid = 0.5*(rbnd+lbnd);
% 
%         while ((rel_err >= rel_err_tol) && (num_gens <= max_gens))
% 
%             gnr_width = mid;
%             gnr_pitch = gnr_width/width_fraction;
%             gnr_space = gnr_pitch - gnr_width;
% 
%             [delay_top delay_side R_top R_top_alt R_side L_vec C_gnr_vec C_gnr_raw_vec Nch_vec mfp_eff_vec ] = ...
%             xcm.calc_gnr_params_combined_multiple_widths( ...
%                 num_layers, gnr_width, gnr_space, gnr_length, temp_K, mfp_defect, rho_interlayer, prob_backscattering, ...
%                 Ef, contact_resistance, epsrd, height_dielectric );
% 
% 
%             delay_gnr = delay_top;
%             pass = (delay_gnr <= delay_max);
% 
%             if(pass)
%                 last_valid_width = mid; % save current value
%                 rbnd = mid; % decrease wire width and continue searching
% 
%                 best_R = R_top;
%                 best_C = C_gnr_vec;
%                 best_C_pul = C_gnr_vec/gnr_length;
%                 best_delay = delay_gnr;
% 
%                 % Only update the error when we have a passing solution, so we
%                 % don't accidentally accept a failing solution
%                 err = delay_max - delay_gnr;
%                 rel_err = abs(err)/delay_max;
% 
%             else % fail
%                 if (last_valid_width == -1)
%                     % We've never found a valid width - need to expand search to higher regimes (rbnd not guaranteed to pass)
%                     lbnd = mid;
%                     rbnd = 2*rbnd;
%                 else % just home in on the line between passing and failing (rbnd gauranteed to pass)
%                     lbnd = mid;
%                 end
%             end
% 
%             mid = 0.5*(lbnd+rbnd); % reset midpoint and keep on testing
%             num_gens = num_gens + 1;
% 
%         end
%     end
% end
% 
% gnr_width = last_valid_width;
% gnr_pitch = gnr_width/width_fraction;
% gnr_delay = best_delay;
% R = best_R;
% C = best_C;
% C_pul = best_C_pul;

            