function [L T N gfrac_L gfrac_T] = find_LT_combination(L,T,slack)
% Finds a suitable combination of L and T to give us an integer number of
% unit cells per side
% L - nominal chip length (in gate pitches)
% T - nominal cell period (in gate pitches)
% slack - fraction (0-1) that tells us when to try increasing/decreasing L,
%   and when to start increasing the cell period instead
%   If L mod T is within slack of a unit cell, we just cut or add some
%   gates.
%   If L mod T is NOT within slack of a unit cell, we increase the cell
%   period until we get closer.

r = mod(L,T);
Lorig = L;
Torig = T;

while r > 0
    if (r <=  slack*T)
        L = L - r;
    elseif (r >= (1-slack)*T)
        L = L + (T-r);
    else
        T = T+1;
    end
    r = mod(L,T);
end

N = L/T;
gfrac_L = (L-Lorig)/Lorig;
gfrac_T = (T-Torig)/Torig;
