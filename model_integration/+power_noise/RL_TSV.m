%%%%%% TSV resistance, inductance and capacitance calculation
%%%%%% Li Zheng, 9/21/2012

function [R,L] = RL_TSV(H,D,barrier_thickness,rho_barrier,rho,mu,TSVpitch,Ppitch,temperature_K)

%%%%%% Test Codes %%%%%
% clear
% clc
% close all
% rho = 1.68*10^-8;                                 %copper resistivity
% mu = 1.257e-6;                                    %copper permeability
% H = 50e-6;
% D = 7e-6;
% TSVpitch = 300e-6;  %TSV pitch (P and G)
% Ppitch = 424e-6;    %power or ground pad pitch
%%%%%%%%%%%%%%%%%%%%%%%%%%

resistivity_bulk = 17.2e-9;
electron_mfp = 39e-9; % (m) Mean free path of electrons in copper
specularity_coeff = 0.55;
reflection_coeff = 0.43;
width = D;
height = D;
wire_length = H;
%[R rho_cu] = xcm.calc_cu_wire_resistance_size_dependent(resistivity_bulk,width,height,wire_length,electron_mfp,specularity_coeff,reflection_coeff); % Use size-dependent resistivity, rather than bulk
[R, rho_cu, R_cu, R_barrier] = xcm.calc_cu_wire_resistance_size_dependent(resistivity_bulk,width,height,barrier_thickness,rho_barrier,wire_length,electron_mfp,specularity_coeff,reflection_coeff,temperature_K);

R = rho_cu*H/(pi*(D/2)^2);
Ls = mu*H/(2*pi)*log(1+2.84*H/(pi*(D/2)));
Lm = 4*0.199*mu*H*log(1+1.0438*H/TSVpitch)-4*0.199*mu*H*log(1+1.0438*H/Ppitch);
L = Ls+Lm;
end
