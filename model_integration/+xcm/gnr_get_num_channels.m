function [ Nch] = gnr_get_num_channels( Ef, gnr_width, T )
% Calculates number of conducting channels in a GNR
% Based on code from Vachan Kumar

%ld = 1000;
% % Ef = 0.3;
q = 1.6e-19; %% electron charge in coulomb
h=6.62e-34; %%% Plancks constant in J-s
vf=8e5; %% Fermi velocity in m/s
kT=0.0258*q*T/300; %% Thermal energy in Joules
%%% calculate the number of conduction channels 
%%% Only Ef+20kT is considered as the upper bound
deltaE=2./(gnr_width*1e9); %% width is in m
Nmax_pos = floor((Ef+20*0.0258*T/300)./deltaE-1/3);
Nmax_neg = floor((Ef+20*0.0258*T/300)./deltaE+1/3);

if (Nmax_pos < 0)
    Nmax_pos = 0; 
    Narr = -(Nmax_neg):1:Nmax_pos;
    Ntot=Nmax_pos+abs(Nmax_neg)+1;
else
Narr = -(Nmax_neg):1:(Nmax_pos);
Ntot=Nmax_pos+abs(Nmax_neg)+1;
end

%% Discretized method -- slightly faster, but not faster **enough** to switch over.

Esub = h*vf./(2*gnr_width) * abs(Narr+1/3); %% sub-band energy in Joules
Nch_pre = 1./(1+exp((Esub-Ef*q)./kT) );

%% outputs
Nch = sum(Nch_pre);


