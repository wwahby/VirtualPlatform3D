function [R rho_cu] = calc_cu_wire_resistance_size_dependent(resistivity_bulk,width,height,length,electron_mfp,specularity_coeff,reflection_coeff)

[rho_cu] = cu_resistivity_fsms(resistivity_bulk,film_thickness,grain_size,electron_mfp,specularity_coeff,reflection_coeff);

R = rho_cu*length/width/height;