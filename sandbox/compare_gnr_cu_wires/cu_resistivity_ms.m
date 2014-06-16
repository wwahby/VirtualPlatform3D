function [rho_ms delta_rho_ms] = cu_resistivity_ms(resistivity_bulk,grain_size,electron_mfp,reflection_coeff)
% Calculates resistivity increase in copper due to surface scattering
% Uses Mayadas and Shatzkes model as described in Tik Sun PhD Thesis 2009
% Typically in small wires grain size is on the order of the thinnest wire
% dimension

R = reflection_coeff;
alpha = electron_mfp/grain_size * R/(1-R);
mod_factor = 3*(1/3-1/2*alpha + alpha^2 - (alpha^3)*log(1+1/alpha));

rho_ms = resistivity_bulk / mod_factor;

delta_rho_ms = rho_ms - resistivity_bulk;