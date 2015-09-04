function [delay_rc, R_gnr, C_gnr] = gnr_calc_delay_rc_kumar(R_source, C_source, C_load, R_contact, epsr_dielectric, space_vertical, gnr_width, gnr_length, mfp_eff, prob_backscattering, num_layers, rho_interlayer, Ef, temp_K)


Rq = 12.9e3; % (Ohm) Quantum resistance of one GNR layer


% Number of conductive channels
Nch = xcm.gnr_get_num_channels( Ef, gnr_width, temp_K );

% Parasitic resistance at each end of the GNR XC
Rp = 0; % Kumar method includes contact and quantum resistance already.

% GNR Resistance
contact_resistance = R_contact;
N_segs_kumar = 10;
[Rnet, Reff] = xcm.calc_gnr_resistance_kumar(num_layers, gnr_width, gnr_length, temp_K, mfp_eff, rho_interlayer, prob_backscattering, Ef, contact_resistance, N_segs_kumar);



% GNR Capacitance per unit length
[C_top, C_side, Ctot_lumped, Cc_lumped, Ctot_lumped_no_coupling] = xcm.gnr_get_capacitance_ns2015(epsr_dielectric, space_vertical, space_vertical, num_layers, Ef, gnr_width, temp_K);

R_gnr = Rnet; %R_tc; % top-contacted GNR resistance (ignoring parasitic resistances at interfaces)
C_gnr = C_top*gnr_length;

delay_rc = xcm.calc_delay_elmore( R_source, C_source, C_load, Rp, R_gnr, C_gnr);
%delay_rc = 0.69*( R_source*(C_source + C_load) + Rp*C_load + R_gnr*C_load + (R_source + Rp/2)*C_gnr) + 0.38*R_gnr*C_gnr;