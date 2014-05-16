function [use_gnr gnr_width gnr_pitch gnr_delay] = find_best_gnr_interconnect(num_layers, gnr_length, delay_cu, width_fraction, pitch_orig, temp_K, mfp_defect, rho_interlayer, prob_backscattering, Ef, contact_resistance, epsrd, height_dielectric)

max_gens = 8; % maximum number of generations for the binary search


gnr_width = pitch_orig*width_fraction;
gnr_space = pitch_orig - gnr_width;

% First check that the GNR delay is better than Cu at copper dimensions
[delay_top delay_side R_top R_top_alt R_side L_vec C_gnr_vec C_gnr_raw_vec Nch_vec mfp_eff_vec ] = ...
    xcm.calc_gnr_params_combined_multiple_widths( ...
        num_layers, gnr_width, gnr_space, gnr_length, temp_K, mfp_defect, rho_interlayer, prob_backscattering, ...
        Ef, contact_resistance, epsrd, height_dielectric );


if (delay_top > delay_cu) % GNR no better than Cu -- don't bother trying to shrink the wires
    use_gnr = 0;
    num = 1;
    denom = 1;
    pass = 0;
    fail = 1;
else
    num = 1;
    denom = 1;
    pass = 1;
    fail = 0;
    pitch_best = pitch_orig;
    use_gnr = 1;
    
    %% Binary search to find best dimensions we can get away with
    for gen = 1:max_gens

        % update fraction of original pitch to find new pitch to try
        num = num*2 - pass + fail;
        denom = 2*denom;
        
        pitch_new = num/denom * pitch_orig;
        gnr_width = pitch_new * width_fraction;
        gnr_space = pitch_new - gnr_width;
        
        
        [delay_top delay_side R_top R_top_alt R_side_vec L_vec C_gnr_vec C_gnr_raw_vec Nch_vec mfp_eff_vec ] = ...
            xcm.calc_gnr_params_combined_multiple_widths( ...
                num_layers, gnr_width, gnr_space, gnr_length, temp_K, mfp_defect, rho_interlayer, prob_backscattering, ...
                Ef, contact_resistance, epsrd, height_dielectric );
            
        % Did we pass or fail?
        pass = 1*(delay_top <= delay_cu);
        fail = 1*(delay_top > delay_cu);
        
        if (pass == 1)
            pitch_best = pitch_new;
            delay_best = delay_top;
            width_best = gnr_width;
            
        end
    end
    
end

if (use_gnr == 1)
    gnr_pitch = pitch_best;
    gnr_delay = delay_best;
    gnr_width = width_best;
else
    gnr_width = 0;
    gnr_delay = 0;
    gnr_pitch = 0;
end

            
        
            
            
            
            
            