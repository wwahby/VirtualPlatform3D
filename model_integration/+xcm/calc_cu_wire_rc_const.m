function [tau_rc_cu R_wire C_wire rho_cu C_wire_pul] = calc_cu_wire_rc_const(resistivity_bulk,width,height,wire_length,electron_mfp,specularity_coeff,reflection_coeff,epsr_dielectric,horiz_space,vert_space)

[R_wire rho_cu] = xcm.calc_cu_wire_resistance_size_dependent(resistivity_bulk,width,height,wire_length,electron_mfp,specularity_coeff,reflection_coeff);
[C_wire C_wire_pul] = xcm.calc_cu_wire_capacitance(epsr_dielectric,wire_length,width,height,horiz_space,vert_space);

tau_rc_cu = R_wire * C_wire;