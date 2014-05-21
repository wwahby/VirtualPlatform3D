function [R rho_cu] = calc_cu_wire_resistance_size_dependent(resistivity_bulk,width,height,wire_length,electron_mfp,specularity_coeff,reflection_coeff)

film_thickness = min(width,height);
grain_size = film_thickness;
[rho_cu] = xcm.cu_resistivity_fsms(resistivity_bulk,film_thickness,grain_size,electron_mfp,specularity_coeff,reflection_coeff);

R = rho_cu*wire_length/width/height;