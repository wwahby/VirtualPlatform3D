function [C_wire C_wire_pul] = calc_cu_wire_capacitance(epsr_dielectric,wire_length,width,height,horiz_space,vert_space)

[cap_const cg cm] = xcm.calc_capacitance_constant_full(width,horiz_space,height,vert_space);

eps0 = 8.854e-12; % (F/m) Vacuum permittivity
epsd = eps0*epsr_dielectric;

C_wire_pul = epsd * cap_const;
C_wire = C_wire_pul * wire_length;