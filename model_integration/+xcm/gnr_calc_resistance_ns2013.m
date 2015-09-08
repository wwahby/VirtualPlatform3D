function [ R_tc, R_tc_pul, R_sc, R_sc_pul ] = gnr_calc_resistance_ns2013(gnr_width, gnr_length, mfp_eff, num_layers, rho_interlayer, Ef, temp_K, N_segs)
% Calculates resistance of a MLGNR interconnect
% == INPUTS ==
% gnr_width (m) : Width of the graphene nanoribbon
% gnr_length (m) : length of the interconnect
% mfp_eff (m) : Effective mean free path in the interconnect. Typical is 300nm, optimistic is 1000nm.
% num_layers (-) : Number of graphene layers in the MLGNR XC
% rho_interlayer (Ohm*m) : resistivity in the vertical direction of the MLGNRXC. Typical values range from 0.3-30 Ohm*cm
% Ef (eV) : Fermi level in the MLGNR. 0 for undoped, ~0.2 for doped.
% temp_K (K) : Ambient temperature
% N_segs (-) : Number of segments to split the interconnect into. Most accurate for N_segs = 4.

%% Constants
delta = 0.35e-9; % (m) Interlayer spacing in MLGNRs
Rq = 12.9e3; % (Ohm) Quantum resistance of one GNR layer
dx = gnr_length/N_segs;

%% Important quantities
Nch = xcm.gnr_get_num_channels( Ef, gnr_width, temp_K );

Ry = rho_interlayer*delta/gnr_width/dx; % (Ohm)
Rx = Rq/Nch/mfp_eff * dx; % (Ohm)

%% Side-contacted resistance
R_sc_seg = Rx/num_layers;

%% Determine resistance
R_tc_seg = Rx;
for nn = 2:num_layers
    R_tc_seg = 1/(1/(R_tc_seg + Ry) + 1/(Rx));
end

R_sc = R_sc_seg * N_segs;
R_tc = R_tc_seg * N_segs;

R_sc_pul = R_sc/gnr_length;
R_tc_pul = R_tc/gnr_length;