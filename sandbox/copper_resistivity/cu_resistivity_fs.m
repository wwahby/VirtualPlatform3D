function [rho_fs delta_rho_fs] = cu_resistivity_fs(resistivity_bulk,film_thickness,electron_mfp,specularity_coeff)
% Calculates resistivity increase in copper due to grain boundary
% scattering as dimensions shrink
% Uses Fuchs and Sondheimer model as described in Tik Sun 2009 PhD Thesis

k = film_thickness/electron_mfp;
p = specularity_coeff;

integrand = @(t) (t.^-3 - t.^-5) .* (1-exp(-k*t))./(1-p*exp(-k*t));

integrated = integral(integrand,1,Inf);

mod_factor = 1/( 1 - (3/2/k)*(1-p)*integrated);

rho_fs = resistivity_bulk*mod_factor;
delta_rho_fs = rho_fs - resistivity_bulk;