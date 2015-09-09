function [width, delay] = find_largest_width_for_wire( delay_func, delay_target, guess_init, delay_tolerance)
% Finds the largest width that a wire can be and still meet a given delay
% target.
%
% == INPUTS ==
% delay_func:       handle to a function that takes wire width (m) as the sole input and
%                       returns wire delay. Delay func should already
%                       incorporate the impact of wire length
% delay_target:     (s) Target for wire delay
% guess_init:       (m) Initial guess for wire width
% delay_tolerance:  (-) Maximum fractional error that can be tolerated in
%                       final delay. E.g. 0.01 = 1% error.
%
% == OUTPUTS ==
% width:    (m) Wire width which satisfies delay criteria
% delay:    (s) Actual wire delay corresponding to width

search_factor = 10;
absolute_min_bound = -inf;
absolute_max_bound = inf;
max_gens_init = 10;
max_gens_bin = 10;
is_increasing = false;

% use this for testing the function
% wire_length = 1000e-6;
% delay_func = @(wire_width) test_delay_func( wire_width, wire_length, wire_width, wire_width);

try
    [width, delay] = misc.positive_magnitude_binary_search( ...
                    delay_func, is_increasing, delay_target, guess_init, delay_tolerance, ...
                    search_factor, absolute_min_bound, absolute_max_bound, ...
                    max_gens_init, max_gens_bin );
catch ME
    % If we couldn't find a proper wire width, throw an exception
    if strcmp(ME.identifier, 'ModifiedSearch:InitialSearch:OutOfBounds')
        % some stuff
        msgid = 'XCM:WireWidthFinder:WidthNotFound';
        msgtxt = 'Wire width not found within search bounds! Try increasing search interval (absolute_min/max_bounds).';
        ME2 = MException(msgid, msgtxt);
        throw(ME2)
    elseif strcmp(ME.identifier, 'ModifiedSearch:InitialSearch:BoundsNotFound')
        % some other stuff
        msgid = 'XCM:WireWidthFinder:WidthNotFound';
        msgtxt = 'Wire width search interval could not be defined during initial search! \nDelay target may be unreachable, or max_gens_init may be set too low. ';
        ME2 = MException(msgid, msgtxt);
        throw(ME2)
    else % If something else went wrong, just rethrow the exception
        rethrow(ME)
    end
end
end


function delay = test_delay_func(xc_width, xc_length, xc_space, space_vertical)
% Test function, not necessary to use this normally

%% ITRS2010 9.5nm node params

Rsd_w = 120; % Ohm/micron
Cg_w = 0.5e-15; % F/micron
Lg = 8.9e-9;
Wmin_p = 3*Lg;
Wmin_n = Lg;

driver_size_mult = 1;
Wp = driver_size_mult * Wmin_p;
Wn = driver_size_mult * Wmin_n;

Rsd_p = Rsd_w * Wp;
Cg_p = Cg_w * Wp;
Cg_n = Cg_w * Wn;

Cg = Cg_p + Cg_n;

%% Driver/Load Parameters
R_source = Rsd_p; % (Ohms)
C_source = Cg; % (F)
C_load = C_source; % (F)
R_contact_gnr = 4.3e3; % (Ohms)
epsr_dielectric = 1.85; % (-)

resistivity_bulk_cu = 17e-9;
mfp_electron_cu = 29e-9;
specularity_coeff = 0.55;
reflection_coeff = 0.45;
horiz_space = xc_space;
R_contact_cu = 100; % random guess for contact resistance of cu to cu
cu_aspect_ratio = 2;
cu_wire_height = cu_aspect_ratio * xc_width;

% Use these for similar values as SB case
% R_source = 8e3;
% C_source = 1e-15*1e6*Wmin_p;
% C_load = 1e-15*1e6*Wmin_p;

%%
[tau_rc_cu, R_wire_cu, C_wire_cu, rho_cu, C_wire_cu_pul] = xcm.calc_cu_wire_rc_const( ...
    resistivity_bulk_cu, xc_width, cu_wire_height, xc_length, ...
    mfp_electron_cu, specularity_coeff, reflection_coeff, ...
    epsr_dielectric, xc_space, space_vertical );

    delay = xcm.calc_delay_elmore( R_source, C_source, C_load, R_contact_cu, R_wire_cu, C_wire_cu);
end
