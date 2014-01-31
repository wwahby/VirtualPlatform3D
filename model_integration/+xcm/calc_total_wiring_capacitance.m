function [Cxc Ltot Cn] = calc_total_wiring_capacitance(pn,Ln,Iidf,l,eps_ild,w_gate)

%% Estimate total wiring capacitance
eps0 = 8.854e-12; % (F/m)

LIDF = l.*Iidf;

pp = [pn(2:end) inf];
pm = [inf pn(1:end-1)];
Cn = 1/2*eps_ild*eps0*(4+pn./pp + pn./pm);

lm = 0;
Ltot = zeros(1,length(Ln));
for lind = 1:length(Ln)
    ln = Ln(lind);
    Ltot(lind) = sum(LIDF(lm+1:ln));
    lm = ln;
end
Cxc = sum(Cn.*Ltot.*w_gate);