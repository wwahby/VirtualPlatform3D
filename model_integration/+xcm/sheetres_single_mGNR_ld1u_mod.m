function [ Nch,Rsheet,leff] = sheetres_single_mGNR_ld1u_mod( Ef,W,P,T,ld)

%ld = 1000;
% % Ef = 0.3;
q = 1.6e-19; %% electron charge in coulomb
h=6.62e-34; %%% Plancks constant in J-s
vf=8e5; %% Fermi velocity in m/s
kT=0.0258*q*T/300; %% Thermal energy in Joules
%%% calculate the number of conduction channels 
%%% Only Ef+20kT is considered as the upper bound
deltaE=2./W; %% width is in nm
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

temp1 = zeros(1,Ntot);
temp2 = zeros(1,Ntot);

%% Original method -- integration
for m=1:Ntot
    Esub(m) = h*vf./(2*W*1e-9).*abs(Narr(m)+1/3); %% sub-band energy in Joules
    
    FD = @(E)((1./(1+exp((E-Ef*1.6e-19)./kT)))); %%% Fermi Dirac integral
    FD_diff = @(E)((1./(1+exp((E-Ef*1.6e-19)./kT))).^2./kT.*exp((E-Ef*1.6e-19)./kT)); %% differential of Fermi Dirac integral
    
    kx_sq = @(E)(2*pi*E./(h*vf)).^2-(Narr(m)+1/3)^2.*pi^2./(W*1e-09)^2;
    kn_sq = (Narr(m)+1/3)^2.*pi^2./(W*1e-09)^2;
    
    Trans = @(E)(P*sqrt(kn_sq)./(W.*sqrt(kx_sq(E)))+1/ld).^-1;%% in nm
    Trans_FD =@(E) Trans(E).*FD_diff(E);
    sigma_1D(m) = 1/12.9e3*quadgk(Trans_FD,Esub(m), Ef*q+20*kT);
    temp1(m) = sigma_1D(m);
 
    Nch_pre(m) = 1./(1+exp((Esub(m)-Ef*q)./kT));
    DOS = @(E)(4*E./(h*vf)./sqrt(E.^2-Esub(m).^2)); %%% DOS in 1D GNR
    DOS_FD = @(E) (DOS(E).*FD(E));
    DOS_FD_diff =@(E)(DOS(E).*FD_diff(E)); %% product of DOS and differential of FD... will be needed for diffusivity
    ns_sub(m) = quadgk(DOS_FD,Esub(m), Ef*q+20*kT);

    ns_sub_diff(m)=quadgk(DOS_FD_diff,Esub(m), Ef*q+20*kT);
end

%% Outputs 
sigma_net = sum(sigma_1D);%%% S-nm
res_pu = 1./sigma_net;
Nch = sum(Nch_pre);
ns_tot = sum(ns_sub);
ns_sub_diff_tot = sum(ns_sub_diff);
D = sigma_net*1e-05./(q^2*ns_sub_diff_tot); %%% cm^2/s
mu = sigma_net*1e-05./(q*ns_tot); %%% cm^2/Vs
leff = 12.9./(res_pu.*Nch)*1e3;%%% nm
Rsheet = res_pu.*W; %% in ohms   

%fprintf('Nch: %.3g \t Rsh: %.3g \t leff: %.3g\n', Nch, Rsheet, leff)


%% Discretized method -- slightly faster, but not faster **enough** to switch over.
% Also sacrifices some accuracy since we have to zero out some poles that
% quadgk would otherwise deal with automatically

% Npts_E = 1e2;
% 
% sigma_1D = zeros(1,Ntot);
% Esub = zeros(1,Ntot);
% Nch_pre = zeros(1,Ntot);
% ns_sub = zeros(1,Ntot);
% ns_sub_diff = zeros(1,Ntot);
% for m=1:Ntot
%     Esub(m) = h*vf./(2*W*1e-9).*abs(Narr(m)+1/3); %% sub-band energy in Joules
% 
%     % Integration bounds
%     Emin = Esub(m);
%     Emax = Ef*q + 20*kT;
%     E = linspace(Emin, Emax, Npts_E);
%     
%     % Ef is in eV, E is in J
%     FD = 1./(1+exp( (E-Ef*q)./kT) );
%     FD_diff = (1./(1+exp((E-Ef*q)./kT))).^2./kT.*exp((E-Ef*q)./kT); %% differential of Fermi Dirac integral
%     
%     kx_sq = (2*pi*E./(h*vf)).^2 - (Narr(m)+1/3)^2.*pi^2./(W*1e-09)^2;
%     kn_sq = (Narr(m)+1/3)^2 .* pi^2./(W*1e-09)^2;
% 
%     Trans = (P*sqrt(kn_sq)./(W.*sqrt(kx_sq))+1/ld).^-1;%% in nm
%     Trans( isnan(Trans) ) = 0;
%     Trans_FD = Trans .* FD_diff;
%     sigma_1D(m) = 1/12.9e3*trapz(E, Trans_FD );
%     temp2(m) = sigma_1D(m);
% 
%     Nch_pre(m) = 1./(1+exp((Esub(m)-Ef*q)./kT));
%     
%     DOS = (4*E./(h*vf)./sqrt(E.^2-Esub(m).^2)); %%% DOS in 1D GNR
%     DOS_FD = DOS.*FD;
%     DOS_FD_diff =DOS.*FD_diff; %% product of DOS and differential of FD... will be needed for diffusivity
%     ns_sub(m) = trapz(E,DOS_FD);
% 
%     ns_sub_diff(m) = trapz(E,DOS_FD_diff);
% end

%% outputs
% sigma_net = sum(sigma_1D);%%% S-nm
% res_pu = 1./sigma_net;
% Nch = sum(Nch_pre);
% ns_tot = sum(ns_sub);
% ns_sub_diff_tot = sum(ns_sub_diff);
% D = sigma_net*1e-05./(q^2*ns_sub_diff_tot); %%% cm^2/s
% mu = sigma_net*1e-05./(q*ns_tot); %%% cm^2/Vs
% leff = 12.9./(res_pu.*Nch)*1e3;%%% nm
% Rsheet = res_pu.*W; %% in ohms   

%fprintf('Nch: %.3g \t Rsh: %.3g \t leff: %.3g\n', Nch, Rsheet, leff)

