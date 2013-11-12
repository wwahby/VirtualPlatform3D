%%%%%% TSV resistance, inductance and capacitance calculation
%%%%%% Li Zheng, 9/21/2012

function [R,L] = RL_TSV(H,D,rho,mu,TSVpitch,Ppitch)

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

R = rho*H/(pi*(D/2)^2);
Ls = mu*H/(2*pi)*log(1+2.84*H/(pi*(D/2)));
Lm = 4*0.199*mu*H*log(1+1.0438*H/TSVpitch)-4*0.199*mu*H*log(1+1.0438*H/Ppitch);
L = Ls+Lm;
end
