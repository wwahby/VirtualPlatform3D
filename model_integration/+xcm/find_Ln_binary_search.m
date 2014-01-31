function [Ln Nm] = find_Ln_binary_search(Lm,Nm,lmax,LIDF,LHSF,RHSF,force_pn)
% Lm:   (gate pitches) length of last interconnect routed
% Nm:   (-) Number of XCs of length Lm remaining to be routed
% lmax: (gate pitches) Maximum interconnect length
% LIDF: 1xlmax+1 vector containing L * Iidf(L)
% LHSF: @(LIDF,Lm,Ln) Function handle for the left side of the equation,
%           essentially the area available multiplied by some constants
% RHS:  @(LIDF,Lm,Ln) Function handle for the RHS of the equation
%           essentially the total length or area of the interconnects to be
%           routed
% force_pn: 0 = pn and Ln determined simultaneously. 1 = force pn to be one
%           gate pitch. This is used when we determine that the smallest pn
%           which simultaneously satifies the timing and areal constraints
%           is actually smaller than the gate pitch. Generally the gate
%           pitch should be the minimum achievable, so we
%           can't have our interconnects being placed on a smaller pitch.


% First, let's check to see that it's possible to find ANY Ln such that an
%   entire set of interconnects of a particular length is routed
% The smallest Ln that could be routed will obviously be Ln = Lm+1, so
% we'll see if that would work
LHS = LHSF(LIDF,Lm,Lm+1,Nm);
RHS = RHSF(LIDF,Lm,Lm+1,Nm);
if (RHS > LHS) % Can't route the Nm remaining Lm's AND all the XCs of length Ln
    
    % Can we route all of the remaining Lm's at least?
    LHS = LHSF(LIDF,Lm,Lm,Nm);
    RHS = RHSF(LIDF,Lm,Lm,Nm);
    
    if(RHS > LHS) % So we can't even route all the Lm's!
        
        % Figure out how many XCs we CAN route, and call it a day
        Ln = Lm; % Last XC routed is STILL length Lm
        
        % The form of this depends on whether we're forcing pn=1 or not
        if (force_pn == 0)
            Nr = floor(LHS/Ln/Lm); % routable XCs
        else
            Nr = floor(LHS/Lm);
        end
        
        Nm = Nm - Nr; % Remaining XCs of length Lm
        
    else % We CAN route all the Lm's, just not all the Ln's
        
        % Figure out how many Ln XCs we can route
        Ln = Lm+1; % We already found that we couldn't route ALL of the Ln's, but we can route all of the Lm's
        
        % The form of this depends on whether we're forcing pn=1 or not
        if (force_pn == 0)
            A_rem = LHS - Ln*Lm*Nm;
            Nr = floor(A_rem/Ln^2);
        else
            A_rem = LHS - Lm*Nm;
            Nr = floor(A_rem/Ln);
        end
         
        Nm = LIDF(Ln+1)/Ln - Nr; % Ln+1 because of indexing length(1) = 0 in our case, and divide LIDF by Ln to get just Iidf.
    end
else
    % If we get past all those edge cases, we should actually be able to
    % find an Ln for which all XCs of length Ln (and any remaining Lm's)
    % are routed


    %% Binary search for Ln
    rbnd = lmax;
    lbnd = Lm+1;
    % if (Nm == 0) % no remaining XCs to route in layer m, start with m+1
    %     lbnd = Lm+1;
    % else % finish routing XCs from layer m
    %     lbnd = Lm;
    % end

    found_Ln = 0;
    
    if(rbnd-lbnd <= 1)
        Ln = rbnd;
    else
        while( (rbnd-lbnd > 1) && (found_Ln == 0))
            % Note we're rounding, so if rbnd = lbnd+1, when we calculat Ln we'll
            % get Ln = rbnd. When we get to that situation we'll need to check to
            % see if we can route all of the XCs of length rbnd, or if we can only
            % route some of them
            Ln = round((lbnd + rbnd)/2); % Start at the midpoint.

            % Set RHS so we can compare against LHS
            % Everything is fixed, except for Ln, which changes as lbnd and rbnd
            %   are recalculated
            LHS = LHSF(LIDF,Lm,Ln,Nm);
            RHS = RHSF(LIDF,Lm,Ln,Nm);

            if (LHS == RHS)
                found_Ln = 1;
            elseif (RHS < LHS) % move Ln further right
                if (Ln ~= rbnd)
                    lbnd = Ln;
                else
                    found_Ln = 1; % if we're already at rbnd, then Ln will just be rbnd
                    %Ln = rbnd;
                end
            else % RHS > LHS -- move Ln further left
                if (Ln-lbnd > 1)
                    %lbnd = lbnd;
                    rbnd = Ln;
                else
                    found_Ln = 1;
                    Ln = lbnd;
                end
            end

           Nm = 0;
        end
    end
end

%dstr = sprintf('force_pn: %d\tLm: %d\tLn %d\tNm %d\tlmax: %d\tLHS: %d\tRHS: %d',force_pn,Lm,Ln,Nm,lmax,LHS,RHS);
%disp(dstr)
