function [Cxc Ltot Cn] = calc_total_wiring_capacitance2(pn,Ln,Nm,Iidf,eps_ild,w_gate)
% Inputs
% ===================
% pn (GP, 1xnm) - wire pitch at each wiring tier (nm = # wiring tiers)
% Ln (GP, 1xnm) - longest wire routed in each wiring tier
% Nm (-, 1xnm) - number of interconnects of length Ln remaining to route in
%                   the NEXT wiring tier
% Iidf (-, 1xlmax+1) - wirelength distribution (lmax = maximum wirelength,
%                       plus 1 because it starts with length 0)
% eps_ild (-) - Relative permittivity of interlayer dielectric
% w_gate  (m) - gate pitch
% =
% Outputs
% =====================
% Cxc (F) - Total wiring capacitance
% Ltot (GP) - Total wiring length per wiring tier
% Cn (F/m) - Wire capacitance per length of wires in tier n


%% Estimate total wiring capacitance
eps0 = 8.854e-12; % (F/m)

% Iidf starts with length 0
% Construct length vector
lmax = length(Iidf) - 1;
l = 0:lmax;

LIDF = l.*Iidf;

pp = [pn(2:end) inf];
pm = [inf pn(1:end-1)];
Cn = 1/2*eps_ild*eps0*(4+pn./pp + pn./pm);

lm = 0;
Ltot = zeros(1,length(Ln));
for lind = 1:length(Ln)
    ln = Ln(lind);
    Ltot(lind) = sum(LIDF(lm+1:ln));
    
    % Correct for partially-routed XC lengths
    if (lind > 1)
        Ltot(lind) = Ltot(lind) + Ln(lind-1)*Nm(lind-1);
    end
    Ltot(lind) = Ltot(lind) - Ln(lind)*(Iidf(ln) - Nm(lind));
    
    lm = ln;
end
Cxc = sum(Cn.*Ltot.*w_gate);