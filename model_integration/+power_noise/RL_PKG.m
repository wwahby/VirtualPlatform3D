clear
clc
close all

rho = 1.68*10^-8;                               %copper resistivity
mu = 1.257e-6;                                  %copper permeability
TH = 50e-6;         %TSV height
TD = 10e-6;         %TSV diameter
pitch = 500e-6;
pp = sqrt(2*pitch^2);

[R,L] = power_noise.RL_TSV(TH,TD,rho,mu,pp,pitch);
