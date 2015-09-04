function [delay_rc, R_gnr, C_gnr, Rp] = gnr_calc_delay_rc_NS(R_source, C_source, C_load, R_contact, epsr_dielectric, space_vertical, gnr_width, gnr_length, mfp_eff, num_layers, rho_interlayer, Ef, temp_K)


Rq = 12.9e3; % (Ohm) Quantum resistance of one GNR layer


% Number of conductive channels
Nch = xcm.gnr_get_num_channels( Ef, gnr_width, temp_K );

% Parasitic resistance at each end of the GNR XC
Rp = R_contact + Rq/2/Nch;

% GNR Resistance
Nsegs_NS = 4; % NS2013 resistance calculation is most accurate for 4 segments. Not sure why this should be, but it's much faster than Kumar method
[ R_tc, R_tc_pul, R_sc, R_sc_pul ] = xcm.gnr_calc_resistance_ns2013(gnr_width, gnr_length, mfp_eff, num_layers, rho_interlayer, Ef, temp_K, Nsegs_NS);

% GNR Capacitance per unit length
[C_top, C_side, Ctot_lumped, Cc_lumped, Ctot_lumped_no_coupling] = xcm.gnr_get_capacitance_ns2015(epsr_dielectric, space_vertical, space_vertical, num_layers, Ef, gnr_width, temp_K);

R_gnr = R_tc; % top-contacted GNR resistance (ignoring parasitic resistances at interfaces)
C_gnr = C_top*gnr_length;

delay_rc = xcm.calc_delay_elmore( R_source, C_source, C_load, Rp, R_gnr, C_gnr);
%delay_rc = 0.69*( R_source*(C_source + C_load) + Rp*C_load + R_gnr*C_load + (R_source + Rp/2)*C_gnr) + 0.38*R_gnr*C_gnr;