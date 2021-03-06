function [use_gnr gnr_width gnr_pitch gnr_delay R C C_pul] = find_best_gnr_interconnect_with_repeaters(num_layers, gnr_length, delay_max, repeater_fraction, Ro, Co, min_pitch, width_fraction, pitch_orig, temp_K, mfp_defect, rho_interlayer, prob_backscattering, Ef, contact_resistance, epsrd, height_dielectric)
% initialize binary search parameters
max_gens = 32; % maximum number of binary search generations to run
err = inf;
rel_err = inf;
rel_err_tol = 1e-2;

gnr_width = pitch_orig*width_fraction;
gnr_space = pitch_orig - gnr_width;

min_width = width_fraction*min_pitch;

last_valid_width = -1; % start with invalid width


% Calculate repeater delay constant for suboptimal RI
alpha_rep = 1.44 + 0.53*(repeater_fraction + 1/repeater_fraction); % Repeater delay modification factor when using repeater_fraction < 1


% First check that the GNR delay is better than Cu at copper dimensions
[delay_top delay_side R_top R_top_alt R_side L_vec C_gnr_vec C_gnr_raw_vec Nch_vec mfp_eff_vec ] = ...
    xcm.calc_gnr_params_combined_multiple_widths( ...
        num_layers, gnr_width, gnr_space, gnr_length, temp_K, mfp_defect, rho_interlayer, prob_backscattering, ...
        Ef, contact_resistance, epsrd, height_dielectric );

tau_rc_gnr = R_top * C_gnr_vec;
delay_gnr_rep = alpha_rep*sqrt(Ro*Co*tau_rc_gnr); %Optimal repeatered delay

if (delay_gnr_rep > delay_max) % GNR not good enough to meet spec at Cu dimensions-- don't bother trying to shrink the wires, as they'll only get worse
    use_gnr = 0;
    best_R = R_top;
    best_C = C_gnr_vec;
    best_C_pul = C_gnr_vec/gnr_length;
    best_delay = delay_top;
else
    use_gnr = 1;
    
    % if min pitch works, just end here and move along
    [delay_top delay_side R_top R_top_alt R_side L_vec C_gnr_vec C_gnr_raw_vec Nch_vec mfp_eff_vec ] = ...
    xcm.calc_gnr_params_combined_multiple_widths( ...
        num_layers, min_width, gnr_space, gnr_length, temp_K, mfp_defect, rho_interlayer, prob_backscattering, ...
        Ef, contact_resistance, epsrd, height_dielectric );

    tau_rc_gnr = R_top * C_gnr_vec;
    delay_gnr_rep = alpha_rep*sqrt(Ro*Co*tau_rc_gnr); %Optimal repeatered delay
    
    pass = (delay_gnr_rep <= delay_max);
    
    if(pass)
        last_valid_width = min_width;

        best_R = R_top;
        best_C = C_gnr_vec;
        best_C_pul = C_gnr_vec/gnr_length;
        best_delay = delay_gnr_rep;
    else % check the rest of the widths

        last_valid_width = gnr_width;
        num_gens = 1; % binary search generation counter

        best_R = R_top;
        best_C = C_gnr_vec;
        best_C_pul = C_gnr_vec/gnr_length;
        best_delay = delay_gnr_rep;

        lbnd = min_width;
        rbnd = gnr_width;
        mid = 0.5*(rbnd+lbnd);

        while ((rel_err >= rel_err_tol) && (num_gens <= max_gens))

            gnr_width = mid;
            gnr_pitch = gnr_width/width_fraction;
            gnr_space = gnr_pitch - gnr_width;

            [delay_top delay_side R_top R_top_alt R_side L_vec C_gnr_vec C_gnr_raw_vec Nch_vec mfp_eff_vec ] = ...
            xcm.calc_gnr_params_combined_multiple_widths( ...
                num_layers, gnr_width, gnr_space, gnr_length, temp_K, mfp_defect, rho_interlayer, prob_backscattering, ...
                Ef, contact_resistance, epsrd, height_dielectric );


            %delay_gnr = delay_top;
            tau_rc_gnr = R_top * C_gnr_vec;
            delay_gnr_rep = alpha_rep*sqrt(Ro*Co*tau_rc_gnr); %Optimal repeatered delay
            pass = (delay_gnr_rep <= delay_max);

            if(pass)
                last_valid_width = mid; % save current value
                rbnd = mid; % decrease wire width and continue searching

                best_R = R_top;
                best_C = C_gnr_vec;
                best_C_pul = C_gnr_vec/gnr_length;
                best_delay = delay_gnr_rep;

                % Only update the error when we have a passing solution, so we
                % don't accidentally accept a failing solution
                err = delay_max - delay_gnr_rep;
                rel_err = abs(err)/delay_max;

            else % fail
                if (last_valid_width == -1)
                    % We've never found a valid width - need to expand search to higher regimes (rbnd not guaranteed to pass)
                    lbnd = mid;
                    rbnd = 2*rbnd;
                else % just home in on the line between passing and failing (rbnd gauranteed to pass)
                    lbnd = mid;
                end
            end

            mid = 0.5*(lbnd+rbnd); % reset midpoint and keep on testing
            num_gens = num_gens + 1;

        end
    end
end

gnr_width = last_valid_width;
gnr_pitch = gnr_width/width_fraction;
gnr_delay = best_delay;
R = best_R;
C = best_C;
C_pul = best_C_pul;

if(gnr_pitch >= pitch_orig)
    use_gnr = 0;
end
















%%
% max_gens = 8; % maximum number of generations for the binary search
% 
% % 
% gnr_width = pitch_orig*width_fraction;
% gnr_space = pitch_orig - gnr_width;
% 
% % First check that the GNR delay is better than Cu at copper dimensions
% [delay_top delay_side R_top R_top_alt R_side L_vec C_gnr_vec C_gnr_raw_vec Nch_vec mfp_eff_vec ] = ...
%     xcm.calc_gnr_params_combined_multiple_widths( ...
%         num_layers, gnr_width, gnr_space, gnr_length, temp_K, mfp_defect, rho_interlayer, prob_backscattering, ...
%         Ef, contact_resistance, epsrd, height_dielectric );
% 
% 
% if (delay_top > delay_max) % GNR no better than Cu -- don't bother trying to shrink the wires
%     use_gnr = 0;
%     num = 1;
%     denom = 1;
%     pass = 0;
%     fail = 1;
% else
%     num = 1;
%     denom = 1;
%     pass = 1;
%     fail = 0;
%     pitch_best = pitch_orig;
%     use_gnr = 1;
%     
%     %% Binary search to find best dimensions we can get away with
%     for gen = 1:max_gens
% 
%         % update fraction of original pitch to find new pitch to try
%         num = num*2 - pass + fail;
%         denom = 2*denom;
%         
%         pitch_new = num/denom * pitch_orig;
%         gnr_width = pitch_new * width_fraction;
%         gnr_space = pitch_new - gnr_width;
%         
%         
%         [delay_top delay_side R_top R_top_alt R_side_vec L_vec C_gnr_vec C_gnr_raw_vec Nch_vec mfp_eff_vec ] = ...
%             xcm.calc_gnr_params_combined_multiple_widths( ...
%                 num_layers, gnr_width, gnr_space, gnr_length, temp_K, mfp_defect, rho_interlayer, prob_backscattering, ...
%                 Ef, contact_resistance, epsrd, height_dielectric );
%             
%         % Did we pass or fail?
%         pass = 1*(delay_top <= delay_max);
%         fail = 1*(delay_top > delay_max);
%         
%         if (pass == 1)
%             pitch_best = pitch_new;
%             delay_best = delay_top;
%             width_best = gnr_width;
%             
%         end
%     end
%     
% end
% 
% if (use_gnr == 1)
%     gnr_pitch = pitch_best;
%     gnr_delay = delay_best;
%     gnr_width = width_best;
% else
%     gnr_width = 0;
%     gnr_delay = 0;
%     gnr_pitch = 0;
% end
% 
% R = R_top;
% C = C_gnr_vec;
% C_pul = C_gnr_vec/gnr_length;

            
        
            
            
            
            
            