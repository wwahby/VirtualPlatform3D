function [width, pitch, delay, R, C, C_pul, rho, R_cu, R_barrier] = find_best_cu_wire_with_repeaters(delay_max, repeater_fraction, Ro, Co, width_guess, min_pitch, resistivity_bulk, barrier_thickness, resistivity_barrier, wire_length, width_fraction, aspect_ratio, electron_mfp, specularity_coeff, reflection_coeff, epsr_dielectric, temperature_K)
%Ro, Co are the driver output resistance and input capacitance

% initialize binary search parameters
max_gens = 10; % maximum number of binary search generations to run
err = inf;
rel_err = inf;
rel_err_tol = 1e-3;

lbnd = min_pitch*width_fraction;
rbnd = width_guess*2;
mid = 0.5*(rbnd+lbnd);

last_valid_width = -1; % start with invalid width
num_gens = 1; % binary search generation counter

% Calculate repeater delay constant for suboptimal RI
alpha_rep = 1.44 + 0.53*(repeater_fraction + 1/repeater_fraction); % Repeater delay modification factor when using repeater_fraction < 1

%% First, check the min pitch
width = min_pitch*width_fraction;
height = aspect_ratio * width;
horiz_space = width*(1/width_fraction - 1);
vert_space = height;
    
% [tau_rc_cu, R_wire, C_wire, rho_cu, C_wire_pul] = xcm.calc_cu_wire_rc_const( ...
%     resistivity_bulk, width, height, wire_length, ...
%     electron_mfp, specularity_coeff, reflection_coeff, ...
%     epsr_dielectric, horiz_space, vert_space );

[tau_rc_cu, R_wire, C_wire, rho_cu, C_wire_pul, R_cu, R_barrier] = ...
    xcm.calc_cu_wire_rc_const( ...
        resistivity_bulk, width, height, ...
        barrier_thickness, resistivity_barrier, wire_length, ...
        electron_mfp, specularity_coeff, reflection_coeff, ...
        epsr_dielectric, horiz_space, vert_space, temperature_K );

delay_cu_rep = alpha_rep*sqrt(Ro*Co*tau_rc_cu); %Optimal repeatered delay
pass = (delay_cu_rep <= delay_max);
if(pass)
    last_valid_width = width; % save current value
    best_R = R_wire;
    best_R_cu = R_cu;
    best_R_barrier = R_barrier;
    best_C = C_wire;
    best_rho = rho_cu;
    best_C_pul = C_wire_pul;
    best_delay = delay_cu_rep;
else
%% If min pitch doesn't work, converge on the actual pitch
    while ((rel_err > rel_err_tol) && (num_gens <= max_gens))

        width = mid;
        height = aspect_ratio * width;
        horiz_space = width*(1/width_fraction - 1);
        vert_space = height;

%         [tau_rc_cu, R_wire, C_wire, rho_cu, C_wire_pul] = xcm.calc_cu_wire_rc_const( ...
%             resistivity_bulk, mid, height, wire_length, ...
%             electron_mfp, specularity_coeff, reflection_coeff, ...
%             epsr_dielectric, horiz_space, vert_space, temperature );
        
        [tau_rc_cu, R_wire, C_wire, rho_cu, C_wire_pul, R_cu, R_barrier] = ...
            xcm.calc_cu_wire_rc_const( ...
                resistivity_bulk, width, height, ...
                barrier_thickness, resistivity_barrier, wire_length, ...
                electron_mfp, specularity_coeff, reflection_coeff, ...
                epsr_dielectric, horiz_space, vert_space, temperature_K );


        delay_cu_rep = alpha_rep*sqrt(Ro*Co*tau_rc_cu); %Optimal repeatered delay
        pass = (delay_cu_rep <= delay_max);

        if(pass)
            last_valid_width = mid; % save current value
            rbnd = mid; % decrease wire width and continue searching

            best_R = R_wire;
            best_R_cu = R_cu;
            best_R_barrier = R_barrier;
            best_C = C_wire;
            best_rho = rho_cu;
            best_C_pul = C_wire_pul;
            best_delay = delay_cu_rep;

            % Only update the error when we have a passing solution, so we
            % don't accidentally accept a failing solution
            err = delay_max - delay_cu_rep;
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

width = last_valid_width;
pitch = width/width_fraction;
delay = best_delay;
R = best_R;
R_cu = best_R_cu;
R_barrier = best_R_barrier;
C = best_C;
rho = best_rho;
C_pul = best_C_pul;
