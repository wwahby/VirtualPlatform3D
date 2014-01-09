function Rnand = calc_Rnand(Vdd,Vt,epsr_ox,t_ox,mobility)
% Calculates average drive resistance of 2 input nand gate based on BACPAC
% formula (see http://web.eecs.umich.edu/~dennis/bacpac/models/delay.html)
% This is used to determine the average gate delay
% ===
% Vdd (V)
% Vt (V)
% mobility (cm^2/Vs)
% eps_ox (F/m)
% t_ox (m)

eps0 = 8.854e-12; % (F/m)vacuum permittivity
eps_ox = epsr_ox*eps0;
mobility = mobility/1e4; % Convert from cm^2/Vs to m^2/Vs

Rnand = 0.805*Vdd/(1/2*mobility*eps_ox/t_ox*(Vdd-Vt)^2);