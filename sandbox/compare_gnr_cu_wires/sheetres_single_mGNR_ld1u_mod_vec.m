function [ Nch,Rsheet,leff] = sheetres_single_mGNR_ld1u_mod_vec( Ef,W,P,T,ld)
% Ef: (eV) Fermi energy of the GNR
% W: (nm) Width of the GNR
% P: (0-1) Backscattering probability
% T: (K) Temperature
% ld: (nm) Dopant-induced mean free path

num_widths = length(W);

%% Constants
q = 1.6e-19; % (C) electron charge
h=6.62e-34; % (Js) Planck constant
vf=8e5; % (m/s) Fermi velocity
kT=0.0258*1.6e-19*T/300; % (J) Thermal energy

%% Calculate the number of conduction channels 
% Ef+20kT is considered as the upper bound
deltaE=2./W; %% width is in nm
Nmax_pos = floor((Ef+20*0.0258*T/300)./deltaE-1/3);
Nmax_neg = floor((Ef+20*0.0258*T/300)./deltaE+1/3);

Nmax_pos( Nmax_pos < 0) = 0; % zero out any positive entries that are actually negative

Narr = zeros(1,num_widths);
for wind = 1:num_widths
    Narr(wind) = -(Nmax_neg(wind)):Nmax_pos(wind);
end
Ntot = Nmax_pos + abs(Nmax_neg) + 1;

Ntot_max = max(Ntot);
Esub = zeros(Ntot_max,length(W));
for m=1:Ntot
    Esub(m) = h*vf./(2*W*1e-9).*abs(Narr(m)+1/3); %% sub-band energy in Joules
    
    FD = @(E)((1./(1+exp((E-Ef*1.6e-19)./kT)))); %%% Fermi Dirac integral
    FD_diff = @(E)((1./(1+exp((E-Ef*1.6e-19)./kT))).^2./kT.*exp((E-Ef*1.6e-19)./kT)); %% differential of Fermi Dirac integral
    
    kx_sq = @(E)(2*pi*E./(h*vf)).^2-(Narr(m)+1/3)^2.*pi^2./(W*1e-09)^2;
    kn_sq = (Narr(m)+1/3)^2.*pi^2./(W*1e-09)^2;
    
    Trans = @(E)(P*sqrt(kn_sq)./(W.*sqrt(kx_sq(E)))+1/ld).^-1;%% in nm
    Trans_FD =@(E) Trans(E).*FD_diff(E);
    sigma_1D(m) = 1/12.9e3*quadgk(Trans_FD,Esub(m), Ef*q+20*kT);
 
    Nch_pre(m) = 1./(1+exp((Esub(m)-Ef*q)./kT));
    DOS = @(E)(4*E./(h*vf)./sqrt(E.^2-Esub(m).^2)); %%% DOS in 1D GNR
    DOS_FD = @(E) (DOS(E).*FD(E));
    DOS_FD_diff =@(E)(DOS(E).*FD_diff(E)); %% product of DOS and differential of FD... will be needed for diffusivity
    ns_sub(m) = quadgk(DOS_FD,Esub(m), Ef*q+20*kT);

    ns_sub_diff(m)=quadgk(DOS_FD_diff,Esub(m), Ef*q+20*kT);
end
sigma_net = sum(sigma_1D);%%% S-nm
res_pu = 1./sigma_net;
Nch = sum(Nch_pre);
ns_tot = sum(ns_sub);
ns_sub_diff_tot = sum(ns_sub_diff);
D = sigma_net*1e-05./(q^2*ns_sub_diff_tot); %%% cm^2/s
mu = sigma_net*1e-05./(q*ns_tot); %%% cm^2/Vs
leff = 12.9./(res_pu.*Nch)*1e3;%%% nm
Rsheet = res_pu.*W; %% in ohms   
end

