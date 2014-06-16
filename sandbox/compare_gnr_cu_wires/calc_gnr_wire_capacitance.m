function [C_gnr C_gnr_raw cap_const Cqe Calpha] = calc_gnr_wire_capacitance(num_layers,Ef,temp_K,widths,gnr_length,thickness,height_dielectric,epsrd,prob_backscattering,mfp_defect)
% width: (m) Width of GNR. Pitch is 2*width
eps0 = 8.854e-12; % (F/m) Vacuum permittivity
eps = eps0*epsrd; % Premittivity of interlayer dielectric
h0 = 6.626e-34; % (J*s) Planck constant
q = 1.602e-19; % (C) electron charge 
h = h0/q; % (eV*s) Planck constant
vf = 8e5; % (m/s) Fermi velocity in graphene

%% Precompute all the stuff we need to get the quantum capacitance for all possible widths
num_widths = length(widths);
Nch = zeros(1,num_widths);
Rsheet = zeros(1,num_widths);
leff = zeros(1,num_widths);
cap_const = zeros(1,num_widths);
cg = zeros(1,num_widths);
cm = zeros(1,num_widths);
for wind = 1:length(widths)
    W_nm = widths(wind)*1e9;
    [ Nch_w,Rsheet_w,leff_w] = sheetres_single_mGNR_ld1u_mod( Ef,W_nm,prob_backscattering,temp_K,mfp_defect*1e9);
    
    [cc cg cm] = calc_capacitance_constant_full(widths(wind),widths(wind),thickness,height_dielectric);
    
    cap_const(wind) = cc;
    cg(wind) = cg;
    cm(wind) = cm;
    
    Nch(wind) = Nch_w;
    Rsheet(wind) = Rsheet_w;
    leff(wind) = leff_w;
end
    
Cqe = 4*q^2/h0/vf * num_layers * Nch; % (F/m) Quantum capacitance;

%% Compute individual capacitance components
% Capacitance to wires below and above this one
%Cbot = eps*width*length/2 / height_dielectric; % Turns out cbot and ctop
%are the same
Ctop = eps*widths/2 / height_dielectric;

Cside = eps*thickness./widths; % We're assuming that the width is the same as the pitch

Calpha = Ctop + Cside;
%Cbeta = Cbot + Cside; % ctop and cbot are the same, so calpha and cbeta
%are too

ca = 1/2*cap_const * eps;
%ca = (cg + cm) * eps;
% C_gnr_raw = 4*Cqe.*Calpha./(2*Cqe+Calpha)*gnr_length;
C_gnr_raw = 2*1./(1/2./Cqe + 1./Calpha)*gnr_length;
C_gnr = 2*1./(1/2./Cqe + 1./ca)*gnr_length;

