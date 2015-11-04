function [tau_rc_cu, R_wire, C_wire, rho_cu, C_wire_pul, R_cu, R_barrier] = calc_cu_wire_rc_const(resistivity_bulk, width, height, barrier_thickness, resistivity_barrier, wire_length, electron_mfp, specularity_coeff, reflection_coeff, epsr_dielectric, horiz_space, vert_space, temperature_K)
% This is the top-level function which determines the RC delay for your wire.
% It will also output the R, C, copper resistivity, and C per unit length for the wire
%	Calculates RC delay constant for a copper wire.
%	Inputs:
%		Bulk copper resistivity (Ohm*m)	::	~17e-9 Ohm*m
%		Wire width (m)	::	Get this from published data or ITRS. NOT equal to the smallest patternable dimension!
%		Wire height (m)	::	Typically ~ wire_width * 2ish
%		wire length (m)	::	Depends on type of wire (local vs intermediate vs global)
%		electron mean free path (m)	::	Using 39nm (From Sun's paper I think)
%		specularity coefficient (-)	::	Using 0.55 (from Sun's paper)
%		reflection coefficient (-)	::	Using 0.43 (from Sun's paper)
%		Dielectric relative permittivity (-)	::	~3.0 for low-k dielectric. 2.7 for ULK, 3.9 for SiO2
%		Horizontal space between wires (m)	::	Typically equal to wire width
%		Vertical space between wires (m)	::	Typically similar to the wire height
%	Outputs:
%		RC delay constant (s)
%		Wire Resistance (Ohms)
%		Wire Capacitance (F)
%		Effective Resistivity (Ohm*m)	::	Actual resistivity of the wire material for the given specularity and reflection parameters
%		Wire capacitance per unit length (F/m)	::	Just wire capacitance / length


%[R_wire, rho_cu] = xcm.calc_cu_wire_resistance_size_dependent(resistivity_bulk, width, height, wire_length, electron_mfp, specularity_coeff, reflection_coeff);
[R_wire, rho_cu, R_cu, R_barrier] = xcm.calc_cu_wire_resistance_size_dependent(resistivity_bulk, width, height, barrier_thickness, resistivity_barrier, wire_length, electron_mfp, specularity_coeff, reflection_coeff,temperature_K);


[C_wire, C_wire_pul] = xcm.calc_cu_wire_capacitance(epsr_dielectric, wire_length, width, height, horiz_space, vert_space);

tau_rc_cu = R_wire * C_wire;
