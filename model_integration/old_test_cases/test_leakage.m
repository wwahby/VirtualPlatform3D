
w_trans = 32e-9;
Ioff_per_um = 100e-9;

width_um = 1.5*w_trans*1e6;
w_trans_um = w_trans * 1e6;
Ilk_To = Ioff_per_um*width_um;
To = 300; % (K)
Vdd = 1.25;
Vds = Vdd;
Vth = 0.25; % (V)

kb = 1.381e-23; % (J/K)
T = 390; % (K)
q = 1.602e-19; % (C)
Eg = 1.12; % (eV)

phi_th = kb*T/q;
phi_th0 = kb*To/q;
%Ilk_T = Ilk_To * (T/To)*exp(-q*Eg/2/kb*(1/T-1/To))*(exp(-Vds/phi_th)-1)
Ilk_T = Ilk_To * (T/To)*exp(-q*(Vds-Vth)/2/kb*(1/T-1/To))*(exp(-Vds/phi_th)-1)
%Ilk_T = Ilk_To * (T/To)*exp(-Eg/2/phi_th)*exp(Eg/2/phi_th0)*(exp(-Vds/phi_th)-1);
Ilk_T = abs(Ilk_T);

Plk_per_transistor = Ilk_T*Vdd
Plk_per_transistor_per_um = Ilk_T/width_um

Nt = 1e6;
Plk_per_million_transistors = Plk_per_transistor*Nt