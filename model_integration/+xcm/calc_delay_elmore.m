function tau_rc = calc_delay_elmore( R_source, C_source, C_load, R_parasitic, R_wire, C_wire)
% Returns 50% elmore delay.
% Parasitic resistance is the total additional resistance at either end of
% the interconnect (i.e. is contact resistance at one end, plus any other
% weird effects you want to include)

tau_rc = 0.69*( R_source*(C_source + C_load) + R_parasitic*C_load + R_wire*C_load + (R_source + R_parasitic/2)*C_wire) + 0.38*R_wire*C_wire;