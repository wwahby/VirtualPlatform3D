function [width pitch delay R C C_pul rho] = find_best_cu_wire(delay_max, width_guess, min_pitch, resistivity_bulk, wire_length, width_fraction, aspect_ratio, electron_mfp, specularity_coeff, reflection_coeff, epsr_dielectric)

% initialize binary search parameters
max_gens = 32; % maximum number of binary search generations to run
err = inf;
rel_err = inf;
rel_err_tol = 1e-2;

lbnd = min_pitch*width_fraction;
rbnd = width_guess*2;
mid = 0.5*(rbnd+lbnd);

last_valid_width = -1; % start with invalid width
num_gens = 1; % binary search generation counter

%% First, check the min pitch
width = min_pitch*width_fraction;
height = aspect_ratio * width;
horiz_space = width*(1/width_fraction - 1);
vert_space = height;
    
[tau_rc_cu R_wire C_wire rho_cu C_wire_pul] = xcm.calc_cu_wire_rc_const( ...
    resistivity_bulk, width, height, wire_length, ...
    electron_mfp, specularity_coeff, reflection_coeff, ...
    epsr_dielectric, horiz_space, vert_space );

delay_cu = 1.1*tau_rc_cu; % Yeah this is a bit weird -- following Venkatesan.
pass = (delay_cu <= delay_max);
if(pass)
    last_valid_width = width; % save current value
    best_R = R_wire;
    best_C = C_wire;
    best_rho = rho_cu;
    best_C_pul = C_wire_pul;
    best_delay = delay_cu;
else
%% If min pitch doesn't work, converge on the actual pitch
    while ((rel_err >= rel_err_tol) && (num_gens <= max_gens))

        width = mid;
        height = aspect_ratio * width;
        horiz_space = width*(1/width_fraction - 1);
        vert_space = height;

        [tau_rc_cu R_wire C_wire rho_cu C_wire_pul] = xcm.calc_cu_wire_rc_const( ...
            resistivity_bulk, width, height, wire_length, ...
            electron_mfp, specularity_coeff, reflection_coeff, ...
            epsr_dielectric, horiz_space, vert_space );


        delay_cu = 1.1*tau_rc_cu; % Yeah this is a bit weird -- following Venkatesan.
        pass = (delay_cu <= delay_max);

        if(pass)
            last_valid_width = mid; % save current value
            rbnd = mid; % decrease wire width and continue searching

            best_R = R_wire;
            best_C = C_wire;
            best_rho = rho_cu;
            best_C_pul = C_wire_pul;
            best_delay = delay_cu;

            % Only update the error when we have a passing solution, so we
            % don't accidentally accept a failing solution
            err = delay_max - delay_cu;
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
C = best_C;
rho = best_rho;
C_pul = best_C_pul;
