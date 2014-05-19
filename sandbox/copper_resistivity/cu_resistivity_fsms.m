function [rho_cu delta_rho_fs delta_rho_ms] = cu_resistivity_fsms(resistivity_bulk,film_thickness,grain_size,electron_mfp,specularity_coeff,reflection_coeff)

[rho_fs delta_rho_fs] = cu_resistivity_fs(resistivity_bulk,film_thickness,electron_mfp,specularity_coeff);
[rho_ms delta_rho_ms] = cu_resistivity_ms(resistivity_bulk,grain_size,electron_mfp,reflection_coeff);

% add resistivities as per mathiessen's rule
rho_cu = resistivity_bulk + delta_rho_fs + delta_rho_ms;