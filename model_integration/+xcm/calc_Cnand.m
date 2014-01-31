function Cnand = calc_Cnand(L_gate_min,epsr_ox,t_ox)
% calculate capacitance of minimum-sized nand gate
% L_gate - (m) gate length
% eps_ox - (F/m) gate oxide permittivity]
% t_ox - (m) gate oxide thickness

eps0 = 8.854e-12; % (F/m)vacuum permittivity
eps_ox = epsr_ox*eps0;

Cnand = 2.5*L_gate_min^2*eps_ox/t_ox;