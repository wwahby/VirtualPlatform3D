function [Plk_tot, Plk_per_transistor_per_um] = calc_transistor_leakage_at_temp(temperature_C, leakage_reference_temperature_K, Vdd, Ioff_per_um, w_trans_min_m, transistor_widths_um, transistor_numbers)

width_um = 1.5*w_trans_min_m*1e6;
Ilk_To = Ioff_per_um*width_um;
To = leakage_reference_temperature_K; % (K)
Vds = Vdd;

kb = 1.381e-23; % (J/K)
T = temperature_C + 273.15; % (K)
q = 1.602e-19; % (C)
Eg = 1.12; % (eV)

phi_th = kb*T/q;
Ilk_T = Ilk_To * (T/To)*exp(-q*Eg/2/kb*(1/T-1/To))*(exp(-Vds/phi_th)-1);
%Ilk_T = Ilk_To * (T/To)*exp(-Eg/2/phi_th)*exp(Eg/2/phi_th0)*(exp(-Vds/phi_th)-1);
Ilk_T = abs(Ilk_T);

Plk_per_transistor_per_um = Ilk_T/width_um;
Plk_tot = sum( transistor_numbers .* transistor_widths_um .* Plk_per_transistor_per_um );