function [Ln_vec pn_vec pn_orig_vec Nm_vec] = wire_layer_assignment_alt(Iidf,lmax,Ach,chi,rho_m,epsr_d,Tclk,alpha_t)
% Determines which wires get assigned to which metal tiers
% Each metal tier consists of two layers of metal interconnects (sized
% equally), routed orthogonally. For simplicity we assume interconnect
% height and width are the same, for an aspect ratio of 1.
% Ignoring a bunch of parasitic capacitances for now
% === Inputs ===
%   Iidf:   (-) 1x(lmax+1) vector of the interconnect density function
%   Ach:    (gates^2) Chip surface area
%   chi:    (-) Some sort of point to point correction? Blame Joyner
%   rho_m:  (Ohm*m) resistivity of the interconnect metal
%   epsr_d: (F/m) Relative permittivity of the dielectric separating wires
%   beta:   (-) Fraction of the clock frequency which can be consumed just
%               by propagation delay
%   Tclk:   (s) Clock period
%   alpha_t:    (-) modifier for RC constant. Joyner and Venkatesan use 1.1
%
% === Outputs ===
%   Ln: (gates) Vector containing maximum interconnect length in each
%           interconnect tier
%   pn: (gates) Vector of the interconnect pitch in each interconnect tier
%   Nm: (-) Number of interconnects of length Ln remaining to be routed in
%   the next tier

%% Constants
eps0 = 8.854e-12; % (F/m) permittivity of free space

%% Setup
eps = epsr_d*eps0; % Permittivity of XC dielectric

n_n = 2; % number of metal layers in tier n
e_w = 0.4; % 40% routing efficiency. When we incorporate the via blockage modeling from Sekar's work this will become dynamic

area_available = n_n*e_w*Ach; % Area available for XC routing

% Set up vector of possible XC lengths
%lmax = length(Iidf) - 1;
lengths = 0:lmax;

% Lm is the longest interconnect routed in layer n-1. It starts out at 1,
% since we have to start somewhere, and we want to route the smallest
% interconnects in the lowest metal layers to minimize the additional
% interconnect length
Lm = 0;
Nm = 0; % No interconnects left over to route. This will be updated with an actual number if we can only route *some* of the interconnects
        % of a particular length -- if that happens, we'll carry Nm XCs of
        % that length for routing in the next tier

% Precalculate the argument of the summation
LIDF = lengths.*Iidf;
ln_ind = 1; % Start by finding the first Ln.
Ln = -1;
force_pn = 0;

while (Ln < lmax)
    % We're going to use a binary search to find Ln. To do that we need to set
    % up left and right bounds -- we'll start with the entire space, and the
    % bounds will shrink with each subsequent layer.
    
    
    % For the first layer, wiring delay is a smaller component of total
    % delay (~25%). For longer layers assume that wiring delay is the
    % dominant component of the clock time (~90%)
    if ln_ind == 1
        beta = 0.25;
    else
        beta = 0.9;
    end
    
    % Set up the quantity we want to match with the RHS
    % RHS is the sum of LIDF from Lm+1 to Ln (inclusive)
    % We leave Lm out because all those should have been routed in the
    %   previous tier
    % The indices are off by one because MATLAB indexing starts at 1, but
    %   l(1) = 0
    % We're also including any interconnects of length Lm that *weren't*
    %   routed in the previous tier (there should be Nm of them 
    if(force_pn == 0)
        LHSF = @(LIDF,Lm,Ln,Nm) area_available/(2*chi*sqrt(alpha_t*rho_m*eps/beta/Tclk));
        RHSF = @(LIDF,Lm,Ln,Nm) Ln*Lm*Nm + Ln*sum(LIDF(Lm+2:Ln+1));
    else
        LHSF = @(LIDF,Lm,Ln,Nm) area_available/chi;
        RHSF = @(LIDF,Lm,Ln,Nm) Lm*Nm + sum(LIDF(Lm+2:Ln+1));
    end

    [Ln Nm_new] = find_Ln_binary_search(Lm,Nm,lmax,LIDF,LHSF,RHSF,force_pn);

    Ln_vec(ln_ind) = Ln;
    Nm_vec(ln_ind) = Nm_new;
    if (force_pn == 0)
        pn_vec(ln_ind) = 2*Ln*sqrt(alpha_t*rho_m*eps/beta/Tclk);
        pn_orig_vec(ln_ind) = pn_vec(ln_ind);
    
        if(pn_vec(ln_ind) >= 1)
            ln_ind = ln_ind+1;
            Lm = Ln; % start next iteration with min gate length = this iterations longest gate length
            force_pn = 0;
        else
            force_pn = 1;
            pn_orig_vec(ln_ind) = pn_vec(ln_ind);
            pn_vec(ln_ind) = 1;
            %disp('Recalculating Ln for this layer!')
        end
    else
        % if force_pn = 1, we just reran everything to find Ln if pn=1 gate
        %   pitch
        % In that case, we can just take the Ln we found and move on with our lives
        % Be sure to reset force_pn so we can let it change dynamically for
        %   other layers
        ln_ind = ln_ind+1;
        Lm = Ln;
        force_pn = 0;
    end
end
        


