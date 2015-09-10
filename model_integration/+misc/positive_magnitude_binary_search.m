function [input, output] = positive_magnitude_binary_search( func, is_increasing, target, guess_init, tolerance, search_factor, absolute_min_bound, absolute_max_bound, max_gens_init, max_gens_bin)
% Does a binary search to find a value of input for which func(input) is
% within tolerance of target. Before the binary search an initial quick
% search is done to find the region about which func - target changes
% signs. The search interval starts at guess_init, and is expanded by a
% factor of search_factor until the initial bounds are determined
%
% == INPUTS == 
% func : handle to a function that accepts one input and returns one
%           output. This is the function that will be evaluated to find an input
%           which is appropriately close to the target. Func is assumed to
%           be monotonic with respect to input. Input is assumed to always
%           be NONNEGATIVE.
% is_increasing: boolean. True = func increases monotonically with input.
%                         False = func decreases monotonically with input.
% target: target value. We want to find an input which causes func(input)
%           to be within tolerance of target.
% guess_init: Initial guess for the input
% tolerance: Maximum RELATIVE error allows (i.e. 0.01 = 1% deviation from
%           target is allowed).
% search_factor: Factor by which to increase search bounds in initial
%           search stage. Typically 10.
% absolute_min_bound: Minimum bound beyond which input is invalid. Use -inf
%           if no such bound exists.
% absolute_max_bound: Maximum bound beyond which input is invalid. Use inf
%           if no such bound exists.
% max_gens_init: Maximum number of generations to allow initial search
%           stage to run.
% max_gens_bin: Maximum number of generations to allow binary search stage
%           to run.
%
% == OUTPUTS ==
% input: The input value which gets us close enough to target
% output: the actual result of func(input)



    % Find coarse bounds
    [min_bound, max_bound, input_init, output_init, within_tol] = find_search_bounds( func, is_increasing, target, guess_init, tolerance, search_factor, max_gens_init, absolute_min_bound, absolute_max_bound);
    
    % If the initial search got us close enough, just use those values
    % Otherwise, run the binary search
    if (within_tol)
        input = input_init;
        output = output_init;
    else
        [input, output] = binary_search( func, is_increasing, target, min_bound, max_bound, tolerance, max_gens_bin);
    end
end


function [min_bound, max_bound, input_init, output_init, within_tol] = find_search_bounds( func, is_increasing, target, guess_init, tolerance, search_factor, max_gens, absolute_min_bound, absolute_max_bound)
    %% Find bounds
    % Finds bounds for binary search
    % Start with initial guess g0
    % assumes func is monotonic
    % finds an interval within which func(input) crosses target
    % if func(input_i) < target, g(i+1) = A*g(i)
    % if func(input_i) > target, g(i+1) = g(i)/A
    % Stop when func(input_i) flips relation with target (i.e. is less than when it
    % was previously greater than, etc

    %search_factor = 10; % amount to multiply or divide by at each generation

    keep_going = true;
    sign_flipped = false;
    within_tol = false;
    cur_input = guess_init;

    gen_init = 0;
    dir_sign = ~is_increasing;
    prev_dir_sign = ~is_increasing;
    prev_input = cur_input;
    prefactor = 1;

    while(keep_going)

        cur_output = func(cur_input);
        
        abs_err = abs(target - cur_output);
        rel_err = abs(abs_err/target);
        within_tol = (rel_err <= tolerance);

        prev_dir_sign = dir_sign;
        if (cur_output < target)
            dir_sign = is_increasing;
        elseif (cur_output > target)
            dir_sign = ~is_increasing;
        else % cur_output == target -- just do nothing and the rest of the loop will execute properly
            dir_sign = ~is_increasing;
        end
        
        if (dir_sign)
            prefactor = search_factor;
        else
            prefactor = 1/search_factor;
        end
        
        % do nothing if we somehow exactly match the target
        if (cur_output == target)
            prefactor = 1;
        end

        if ( gen_init > 0) % the second run is the first one where we'll know if we can stop
            if ( dir_sign ~= prev_dir_sign)
                sign_flipped = true;
                min_bound = min(prev_input, cur_input);
                max_bound = max(prev_input, cur_input);
            end
        end


        % Stopping conditions
        input_within_bounds = (cur_input > absolute_min_bound) && (cur_input < absolute_max_bound);
        if (~input_within_bounds) % Throw an exception if this happens so we know what happened.
            msgID = 'ModifiedSearch:InitialSearch:OutOfBounds';
            msgtext = sprintf('Could not find appropriate bounds for search variable within given interval!\nValue: %.3g \t min_bound: %.3g \t max_bound: %.3g', cur_input, absolute_min_bound, absolute_max_bound);
            ME = MException(msgID, msgtext);
            throw(ME)
        end
        
        keep_going = (~sign_flipped) && (~within_tol) && (gen_init < max_gens) && input_within_bounds;

        % Debug print -- can be commented out most of the time
        %fprintf('Initial Search Gen %d: \t Input: %.3g \t Output: %.4d\n\n', gen_init, cur_input, cur_output);
        
        if (keep_going)
            prev_input = cur_input;
            cur_input = prefactor*prev_input;
            gen_init = gen_init + 1;
        elseif (~sign_flipped) % haven't found bounds yet
            msgID = 'ModifiedSearch:InitialSearch:BoundsNotFound';
            msgtext = sprintf('Could not find appropriate bounds for search variable within %d search generations! Try increasing max_gens.', max_gens);
            ME = MException(msgID, msgtext);
            throw(ME);
        end
    end
    
    input_init = cur_input;
    output_init = cur_output;
end


%% Binary search
function [input, output] = binary_search( func, is_increasing, target, min_bound, max_bound, tolerance, max_gens)
    left = min_bound;
    right = max_bound;
    num_gens = 0;
    abs_err = inf;
    rel_err = inf;
    
    while ( (rel_err > tolerance) && (num_gens < max_gens) )
        mid = 1/2*(left+right);
        cur_output = func(mid);
        abs_err = abs(target - cur_output);
        rel_err = abs(abs_err/target);

        if (is_increasing)
            if (cur_output > target)
                right = mid;
            elseif (cur_output < target)
                left = mid;
            else
                left = 1/2*(left + mid);
                right = 1/2*(right + mid);
            end
        else
            if (cur_output > target)
                left = mid;
            elseif (cur_output < target)
                right = mid;
            else
                left = 1/2*(left + mid);
                right = 1/2*(right + mid);
            end
        end

        num_gens = num_gens + 1;
        fprintf('Binsrch Gen %d: \t Input: %.3g \t Output: %.4g\n\n',num_gens, mid, cur_output);
    end
    
    input = mid;
    output = cur_output;
end