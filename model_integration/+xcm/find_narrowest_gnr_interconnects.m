function [width_rc, delay_rc, width_rep, delay_rep, width_cu, delay_cu, width_cu_rep, delay_cu_rep] = find_narrowest_gnr_interconnects(delay_target, delay_tolerance, guess_init, repeater_fraction, R_source, C_source, C_load, R_contact_gnr, epsr_dielectric, space_vertical, xc_length, mfp_eff_gnr, num_layers_gnr, rho_interlayer_gnr, Ef_gnr, temp_K)
    
    delay_func = @(wire_width) calc_gnr_delay(R_source, C_source, C_load, R_contact_gnr, epsr_dielectric, space_vertical, wire_width, xc_length, mfp_eff_gnr, num_layers_gnr, rho_interlayer_gnr, Ef_gnr, temp_K);
    [width_rc, delay_rc] = xcm.find_smallest_width_for_wire( delay_func, delay_target, guess_init, delay_tolerance);
    
    delay_func = @(wire_width) calc_gnr_repeatered_delay(R_source, C_source, C_load, R_contact_gnr, epsr_dielectric, space_vertical, wire_width, xc_length, mfp_eff_gnr, num_layers_gnr, rho_interlayer_gnr, Ef_gnr, temp_K, repeater_fraction);
    [width_rep, delay_rep] = xcm.find_smallest_width_for_wire( delay_func, delay_target, guess_init, delay_tolerance);
    
    resistivity_bulk_cu = 17e-9;
    mfp_electron_cu = 39e-9;
    specularity_coeff = 0.55;
    reflection_coeff = 0.45;
    R_contact_cu = 100; % random guess for contact resistance of cu to cu
    cu_wire_aspect_ratio = 2;
    
    delay_func = @(wire_width) calc_cu_delay(R_source, C_source, C_load, R_contact_cu, resistivity_bulk_cu, wire_width, cu_wire_aspect_ratio*wire_width, xc_length, mfp_electron_cu, specularity_coeff, reflection_coeff, epsr_dielectric, wire_width, space_vertical);
    [width_cu, delay_cu] = xcm.find_smallest_width_for_wire( delay_func, delay_target, guess_init, delay_tolerance);
    
    delay_func = @(wire_width) calc_cu_repeatered_delay(R_source, C_source, C_load, R_contact_cu, resistivity_bulk_cu, wire_width, cu_wire_aspect_ratio*wire_width, xc_length, mfp_electron_cu, specularity_coeff, reflection_coeff, epsr_dielectric, wire_width, space_vertical, repeater_fraction);
    [width_cu_rep, delay_cu_rep] = xcm.find_smallest_width_for_wire( delay_func, delay_target, guess_init, delay_tolerance);
end



function delay_rc = calc_gnr_delay(R_source, C_source, C_load, R_contact_gnr, epsr_dielectric, space_vertical, gnr_width, xc_length, mfp_eff_gnr, num_layers_gnr, rho_interlayer_gnr, Ef_gnr, temp_K)
    [delay_rc, R_gnr, C_gnr, Rp] = xcm.gnr_calc_delay_rc_NS(R_source, C_source, C_load, R_contact_gnr, epsr_dielectric, space_vertical, gnr_width, xc_length, mfp_eff_gnr, num_layers_gnr, rho_interlayer_gnr, Ef_gnr, temp_K);
end



function delay_rep = calc_gnr_repeatered_delay(R_source, C_source, C_load, R_contact_gnr, epsr_dielectric, space_vertical, gnr_width, xc_length, mfp_eff_gnr, num_layers_gnr, rho_interlayer_gnr, Ef_gnr, temp_K, repeater_fraction)
    [delay_rc, R_gnr, C_gnr, Rp] = xcm.gnr_calc_delay_rc_NS(R_source, C_source, C_load, R_contact_gnr, epsr_dielectric, space_vertical, gnr_width, xc_length, mfp_eff_gnr, num_layers_gnr, rho_interlayer_gnr, Ef_gnr, temp_K);
    
    Ro = R_source;
    Co = C_load;
    Rint = R_gnr + Rp;
    Cint = C_gnr;
    
    k = repeater_fraction*sqrt(Rint*Cint/2.3/Ro/Co);
    h = sqrt(Ro*Cint/Rint/Co);
    
    delay_rep = k*( 2.3*Ro/h*(Cint/k + h*Co) + Rint/k*(Cint/k + 2.3*h*Co) );
end

function delay_cu_rc = calc_cu_delay(R_source, C_source, C_load, R_contact_cu, resistivity_bulk_cu, xc_width, cu_wire_height, xc_length, mfp_electron_cu, specularity_coeff, reflection_coeff, epsr_dielectric, xc_space, space_vertical)

	[tau_rc_cu, R_wire_cu, C_wire_cu, rho_cu, C_wire_cu_pul] = xcm.calc_cu_wire_rc_const(resistivity_bulk_cu, xc_width, cu_wire_height, xc_length, mfp_electron_cu, specularity_coeff, reflection_coeff, epsr_dielectric, xc_space, space_vertical);
	delay_cu_rc = xcm.calc_delay_elmore( R_source, C_source, C_load, R_contact_cu, R_wire_cu, C_wire_cu);
    
end

function delay_cu_rep = calc_cu_repeatered_delay(R_source, C_source, C_load, R_contact_cu, resistivity_bulk_cu, xc_width, cu_wire_height, xc_length, mfp_electron_cu, specularity_coeff, reflection_coeff, epsr_dielectric, xc_space, space_vertical, repeater_fraction)

	[tau_rc_cu, R_wire_cu, C_wire_cu, rho_cu, C_wire_cu_pul] = xcm.calc_cu_wire_rc_const(resistivity_bulk_cu, xc_width, cu_wire_height, xc_length, mfp_electron_cu, specularity_coeff, reflection_coeff, epsr_dielectric, xc_space, space_vertical);
  
    Ro = R_source;
    Co = C_load;
    Rint = R_wire_cu + 2*R_contact_cu;
    Cint = C_wire_cu;
    
    k = repeater_fraction*sqrt(Rint*Cint/2.3/Ro/Co);
    h = sqrt(Ro*Cint/Rint/Co);
    
    delay_cu_rep = k*( 2.3*Ro/h*(Cint/k + h*Co) + Rint/k*(Cint/k + 2.3*h*Co) );
    
end