function [Ln_vec pn_vec pn_orig_vec A_wires A_vias] = wla_improved_old(Iidf,gate_pitch,min_pitch,layers_per_tier,routing_efficiency,layer_area,rho_m,epsr_d,alpha_t,Beta,Tclk,Rc)
% Simple wire layer assignment algorithm, including via blockage
% Assign wires to layers to satisfy timing and areal constraints as
% described by Venkatesan (see Joyner thesis)
% Assuming wires have unity aspect ratio (i,e, they are pn/2 wide and pn/2
% tall)
% via blockage is included by counting how many wires remain to be routed,
% and assuming a square (pn/2)^2 via will be required on wiring layer n for
% each wire not yet routed.
%==
% Inputs
% ==
% Iidf - wirelength distribution (number of wires of a particular length
%        Wire lengths should be in integer units of gate pitches
% gate_pitch - (m) Distance between logic gates (NOT the width of the gate
%              on an individual transistor!)
% min_pitch - (m) minimum patternable feature pitch
% layers_per_tier - (-) how many metal layers per wiring tier? Wires are
%                   often routed in horizontal and vertical tiers, so
%                   Joyner and Venkatesan use 2 for this (since they want
%                   the dimensions of wires in the H and V tiers to be the
%                   same)
% routing_efficiency - (-) Fraction [0-1] of layer area that your routing tool
%                        actually uses for wire routing. Usually ~0.4-0.6
% layer_area - (m^2) Area available in each layer for routing wires
% rho_m - (Ohm*m) Resistivity of wires
%         If a 1xn vector is input, the first n entries will be used for
%         the first n wiring tiers, and all tiers above n will use the nth
%         entry.
% epsr_d - (-) relative permittivity of interlayer dielectric
% alpha_t - (-) timing constant related to wiring capacitance and delay
%            Venkatesan gives this as 1.1*6.2
% Beta - (-) Fraction [0-1] of clock period that wire delay in each
%         tier can consume. Handles vectors the same way rho_m does
% Tclk - (s) Clock period
% Rc - (Ohm) Contact resistance between tiers. Use 0 to ignore Rc
%          Like rho_m, you can input a vector to have different Rc's for
%          different tiers

chi = 2/3; % conversion factor -- converts between point-to-point wirelength and total net length
eps0 = 8.854e-12;
eps_d = epsr_d * eps0;

l = 0:length(Iidf)-1;
Iidf = round(Iidf); % rounding this because it's nonsensical to route 0.32 of a wire
LIDF = l.*Iidf;
lmax = find(Iidf,1,'last'); % find last length for which we have more than zero wires


n = 1; % start with first wiring tier
Ln = 1; % length of longest wire routed in this tier (GP)
Lm = 0; % length of longest wire routed in previous tier (GP)

% Have to construct LHSF for each iteration
% Update Ln each time around
% Keep track of the difference between the index of Iidf and actual length
% since IIdf(1) corresponds to length=0
while ((Ln < lmax-1) && (Ln > 0))
    
    % Get the value of these constants for each layer we're looking at
    % Each of these can be input as a vector
    % If it's a vector, grab the nth element
    % If n is larger than the last element, grab the last element
    Beta_n = get_nth_or_last(Beta,n);
    rho_m_n = get_nth_or_last(rho_m,n);
    Rc_n = get_nth_or_last(Rc,n);
    
    pnf = @(Ln) sqrt( 4*alpha_t*rho_m_n*eps_d*(Ln*gate_pitch)^2 / (Beta_n*Tclk - alpha_t*Rc_n*eps_d*(Ln*gate_pitch)) );
    A_wires_n = @(Lm,Ln) chi*pnf(Ln)*gate_pitch*sum(LIDF(Lm+2:Ln+1)); % +2 and +1 in LIDF because we need the indices, not the actual lengths
    A_vias_n = @(Ln) layers_per_tier * pnf(Ln)^2 * sum(Iidf(Ln+2:end)); % +2 in Iidf because we need the index, not the length
    A_req_n = @(Lm,Ln) A_wires_n(Lm,Ln) + A_vias_n(Ln);
    
    LHSF = @(Ln) A_req_n(Lm,Ln);
    RHSF = @(ind) layers_per_tier * routing_efficiency * layer_area; % yes, this is independent of the input.
    
    Ln_ind = binary_search_le(LHSF,RHSF,Lm+1,lmax);
    Ln = Ln_ind - 1;
    
    if(Ln_ind == -1)
        err_str = sprintf('WLA Error! Impossible to route wires in layer %d',n);
        disp(err_str);
        Ln = -1; % use this to exit the loop
    else
        pn = pnf(Ln);
        pn_orig = pn; % store off what we *wanted* the original pitch to be
        
        if (pn < min_pitch) % Have to force p_n to min_pitch and find new Ln
            % Fix pn at min_pitch and reset all the function handles
            pnf = @(Ln) min_pitch;
            A_wires_n = @(Lm,Ln) chi*pnf(Ln)*gate_pitch*sum(LIDF(Lm+2:Ln+1)); % +2 and +1 in LIDF because we need the indices, not the actual lengths
            A_vias_n = @(Ln) layers_per_tier * pnf(Ln)^2 * sum(Iidf(Ln+2:end)); % +2 in Iidf because we need the index, not the length
            A_req_n = @(Lm,Ln) A_wires_n(Lm,Ln) + A_vias_n(Ln);
            LHSF = @(Ln) A_req_n(Lm,Ln);
            
            Ln_ind = binary_search_le(LHSF,RHSF,Lm+1,lmax);
            Ln = Ln_ind - 1;
            if(Ln_ind == -1)
                err_str = sprintf('WLA Error! Impossible to route wires in layer %d',n);
                disp(err_str);
                Ln = -1; % use this to exit the loop
            end
            
            % Now that we've confirmed we can actually route this layer,
            % save off the actual wiring pitch
            pn = min_pitch;
        end
        
        pn_vec(n) = pn;
        pn_orig_vec(n) = pn_orig;
        Ln_vec(n) = Ln;
        A_wires(n) = A_wires_n(Lm,Ln);
        A_vias(n) = A_vias_n(Ln);
        n = n+1;
        Lm = Ln;
    end

end




