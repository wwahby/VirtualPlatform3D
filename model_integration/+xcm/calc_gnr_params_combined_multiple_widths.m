function [delay_top_vec delay_side_vec R_top_vec R_top_alt_vec R_side_vec L_vec C_gnr_vec C_gnr_raw_vec Nch_vec mfp_eff_vec ] = calc_gnr_params_combined_multiple_widths(num_layers,gnr_widths,gnr_spaces,gnr_length,Temp_K,mfp_defect,rho_interlayer,prob_backscattering,Ef,contact_resistance,epsrd,height_dielectric)
% Calculates the effective resistance of a graphene nanoribbon interconnect
% GNR Model - Nishad 2014
% INPUTS
% ==============
% num_layers: (-) Number of layers of grpahene in multi-layer GNR interconnect
% width:    (m) Nanoribbon width (can be a vector)
% length:   (m) Nanoribbon length
% Temp_K:     (K) Temperature
% mfp_defect: (m) Mean free path of electrons in the graphene nanoribbon when considering only defects. 100nm - 1um
% rho_interlayer: (Ohm*m) Vertical conductivity in mGNR
% Ef:       (eV) Fermi level in all levels of GNR (assumes doped GNR)
% contact_resistance: (Ohms) Contact resistance between mGNR XC and metal contacts
% .
% OUTPUTS
% =========
% Reff: (Ohms) Effective resistance of the mGNR interconnect
% Rnet: (Ohms) Total resistance (from pad to pad) of the mGNR interconnect

num_widths = length(gnr_widths); % Number of widths to calculate GNR params for
%num_layers = 8; % number of layers of graphene in GNR

%gnr_length = 20e-6; % (m) length of GNR XC
%gnr_width = 7.8e-9; % (m) width of GNR XC
h_vert = height_dielectric;

%rho_interlayer = 3e-1; %(Ohm*m) c-axis resistivity of GNR - Between 3e-3 and 3e-1 Ohm*m
%Ef = 0.2; % (eV) Fermi energy of the bottom GNR layer
d_m = 0.35e-9; % (m) interlayer separation in GNR

Rq = 12.9e3; % (Ohm) Quantum resistance of a single-layer GNR
Rc = contact_resistance; % (Ohm) Contact resistance between contact pads and GNR XC
vf = 8e5; % (m/s) Fermi velocity in graphene
%mfp_defect = 1000e-9; % (m) MFP from defect scattering in graphene on SiO2 (100-300nm realistic, 1um ideal)
%prob_backscattering = 0.0; % (-) Backscattering probability ( 0 for ballistic, 0.2 realistic)
%epsrd = 16.5; % ILD relative permittivity

kb = 1.381e-23; % Boltzmann constant (J/K)
h0 = 6.626e-34; % (J*s) Planck constant
q = 1.602e-19; % (C) electron charge 
h = h0/q; % (eV*s) Planck constant
eps0 = 8.854e-12; % (F/m) Vacuum permittivity

kbT = kb/q*Temp_K; % (eV) combined Boltzmann temp factor
epsd = epsrd*eps0;


%% Calculate number of conducting channels and effective mean free path

Nch_vec = zeros(1,num_widths);
Rsheet_vec = zeros(1,num_widths);
mfp_eff_vec = zeros(1,num_widths);
C_gnr_vec = zeros(1,num_widths);
C_gnr_raw_vec = zeros(1,num_widths);
gnr_cap_const = zeros(1,num_widths);
Cqe_vec = zeros(1,num_widths);
Calpha_vec = zeros(1,num_widths);
R_top_vec = zeros(1,num_widths);
R_top_alt_vec = zeros(1,num_widths);
R_side_vec = zeros(1,num_widths);
L_vec = zeros(1,num_widths);
delay_top_vec = zeros(1,num_widths);
delay_side_vec = zeros(1,num_widths);

for wind = 1:num_widths
    gnr_width = gnr_widths(wind);
    gnr_space = gnr_spaces(wind);
    gnr_width_nm = gnr_width*1e9;
    mfp_defect_nm = mfp_defect * 1e9;

    [ Nch, Rsheet, leff] = xcm.sheetres_single_mGNR_ld1u_mod( Ef,gnr_width_nm,prob_backscattering,Temp_K,mfp_defect_nm);
    
    lambda_eff = leff*1e-9;
    mfp_eff_vec(wind) = leff;
    Rsheet_vec(wind) = Rsheet;
    Nch_vec(wind) = Nch;


    %% Calculate component resistances
    ry = rho_interlayer*d_m/gnr_width/gnr_length;
    rx = Rq/Nch/lambda_eff*gnr_length;

    %% Top contacted MLGNR resistance

    rn = zeros(1,num_layers);
    rn(1) = rx;
    for n = 2:num_layers
        ra = rn(n-1) + ry;
        rb = rx;
        rn(n) = (ra^-1 + rb^-1)^-1;
    end

    r_tc = rn(end);

    %% Side-contacted MLGNR resistance

    r_sc = rx/num_layers;

    %% Total resistance

    R_top = r_tc ;
    R_side = r_sc ;

    [Rnet Reff] = xcm.calc_gnr_resistance_kumar(num_layers,gnr_width,gnr_length,Temp_K,mfp_defect,rho_interlayer,prob_backscattering,Ef,contact_resistance);

    R_top = Reff;
    R_top_alt_vec(wind) = r_tc;
    R_top_vec(wind) = R_top;
    R_side_vec(wind) = R_side;
    %% Inductance

    Lk = h0/4/q^2/vf/num_layers/Nch; %(H/m) Kinetic inductance
    l_pul = Lk;
    %L = Lk;

    %% Capacitance

%     M1 = @(a) 2*pi./( log( 2*(1+(1-a.^2).^(1/4))./(1-(1-a.^2).^(1/4)) ));
%     M2 = @(a) 2/pi*log( (2*(1+sqrt(a)))./(1-sqrt(a)) );
% 
%     M = @(a) ((0 <= a) & (a<= 1/sqrt(2))) .* M1(a) + ...
%              ((1/sqrt(2) < a) & (a<= 1)) .* M2(a);
% 
%     Cq = 4*q^2/h0/vf * num_layers * Nch; % (F/m) Quantum capacitance;
%     Ce = epsd*M(tanh(pi*gnr_width/4/h_vert));
% 
%     c_pul = (Cq^-1 + Ce^-1)^-1;
    
    thickness = 0.5e-9*num_layers;
    [C_gnr C_gnr_raw cap_const Cqe Calpha] = xcm.calc_gnr_wire_capacitance(num_layers,Ef,Temp_K,gnr_width,gnr_space,gnr_length,thickness,height_dielectric,epsrd,prob_backscattering,mfp_defect);
    
%     C_gnr_vec(wind) = C_gnr;
    C_gnr_vec(wind) = C_gnr;
    C_gnr_raw_vec(wind) = C_gnr_raw;
    gnr_cap_const(wind) = cap_const;
    Cqe_vec(wind) = Cqe;
    Calpha_vec(wind) = Calpha;

    %% Delay

    % Get resistances in per-unit-length form
    rtc_pul = Reff/gnr_length;
    rsc_pul = R_side/gnr_length;

    % Get full capacitance and inductance
    C = C_gnr_vec(wind);
    c_pul = C/gnr_length;
    L = l_pul * gnr_length;
    L_vec(wind) = L;

    wn = 1/gnr_length*sqrt(2/l_pul/c_pul);
    alpha_tc = rtc_pul/2/wn/l_pul;
    alpha_sc = rsc_pul/2/wn/l_pul;

    delay = @(wn,gamma) (1.047*exp(gamma/0.85) + 1.39*gamma)/wn;
    delay_side = delay(wn,alpha_sc);
    delay_top = delay(wn,alpha_tc);

    delay_top_kumar = 0.38*Rq*c_pul*gnr_length/Nch + 0.38*rtc_pul*c_pul*gnr_length^2;
    delay_side_kumar = 0.38*Rq*c_pul*gnr_length/Nch + 0.38*rsc_pul*c_pul*gnr_length^2;

    delay_top = delay_top_kumar;
    delay_side = delay_side_kumar;
    
    delay_top = (Rq + rtc_pul*gnr_length).*(c_pul*gnr_length); % RC
    
    delay_top_vec(wind) = delay_top;
    delay_side_vec(wind) = delay_side;
end


