function cap_const = calc_capacitance_constant(ar_wire, width_fraction)
% cap_const = calc_capacitance_constant(ar_wire, width_fraction)
% ==============
% Calculates the multiplicative constant for the capacitance term of the
% interconnect delay
% == INPUTS ==
%   ar_wire: (-) Wire aspect ratio (height/width)
%   width_fraction: (-) wire width as a fraction of wire pitch (width/pitch)
% == OUTPUTS ==
%   cap_const: (-) A multiplicative term to find the effective wire
%                  capacitance (C = cap_const * eps*A/d).
% The reason we have this
% multiplicative factor is because of the wires surrounding the wire of
% interest, and the fact that they may be at Vdd or ground with some
% probability. For more details, see Venkatesan's thesis (eq 2.4-2.7)
% ===
% For 1:1 AR use ar_wire = 1, width_fraction = 0.5;
% For more realistic dimensions for ICs from 2005-2013ish, use ar_wire =
% 1.8, gamma_w = 0.25

%% Helper functions
% All of these variables are given in units of wire pitches
% (they're actually just unitless, since they all appear in ratios, but
% it's easiest to think of them as fractions of the pitch)
% W: wire width
% T: wire height
% S: space between wire and next wire
% H: space between top of wire and next metal layer

cgn = @(W,S,T,H) W/H+1.086*(1+0.685*exp(-T/1.343/S)-0.9964*exp(-S/1.421/H))*(S/(S+2*H))^0.0476*(T/H)^0.337;
cmn = @(W,S,T,H) T/S*(1-1.897*exp(-H/0.31/S)*exp(-T/2.474/S) + 1.302*exp(-H/0.082/S) - 0.1292*exp(-T/1.326/S) ) ...
        + 1.722*(1-0.6548*exp(-W/0.3477/H))*exp(-S/0.651/H);
    
cn = @(W,S,T,H) 2*(cgn(W,S,T,H) + cmn(W,S,T,H));

%% other AR


W = width_fraction;
S = 1-W;
T = ar_wire*W;
H = S;

cap_const = cn(W,S,T,H);