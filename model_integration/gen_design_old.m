function [iidf l Ln pn pn_orig Cxc Ltot Cn Pdyn Plk Pw Prep Ng_act N_tsvs iidf_rewire] = gen_design_old(Ng,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring)

%% Presize the chip and TSVs
Ns = Ng/S;
Lx = round(sqrt(Ns));
Ach_tier_gp = Ach_m2/gate_pitch^2/S; % Force chip to use specified area 
Ach_tier_m2 = Ach_m2/S;

% Recalculate these to make sure everything is a nice integer
Ns = Lx^2;
Ng = Ns*S;

% Size the TSVs
h_tsv = ceil(h_tsv_m/gate_pitch);
w_tsv = ceil(h_tsv/AR_tsv);

%% Size the chip so we have a nicely divisible number of unit cells per side
if ((S == 1) || (use_joyner == 1)) % no TSVs for single layer device, or when using original Joyner model
    w_tsv = 0;
    Lxc = Lx;
    Nsc = Ns;
    Ngc = Ng;
    Nuc_1d = 1;
    N_tsvs = 0;
else
    %Nsp = floor( (1+Atf_max)*Ns );
    Nsp = floor(Ns/(1-Atf_max));
    Lxp = floor(sqrt(Nsp));
    Nsp = Lxp^2;
    Tp = ceil(w_tsv/sqrt(Atf_max));

    slack = 0.2;
    [Lxc Tc Nuc_1d gfrac_L gfrac_T] = find_LT_combination(Lxp,Tp,slack);
    N_tsvs = Nuc_1d^2;
end

Nsc = Lxc^2;
Ngc = Nsc*S;

g_tsv = (Nuc_1d*w_tsv)^2; % number of gates displaced by TSVs
Atf_act = g_tsv/Nsc;

Ns_act = Nsc - g_tsv;
Ng_act = Ns_act*S;

repstr1 = sprintf('Ng_nom: %.4g\tNg_cor: %.4g\tNg_act: %.4g\tAtf_act: %.4g',Ng,Ngc,Ng_act,Atf_act);
disp(repstr1)

%% Calculate WLD
iidf = calc_Iidf_corrected(alpha,k,p,Lx,S,h_tsv,Nuc_1d,w_tsv);
%iidf = calc_Iidf(alpha,k,p,round(sqrt(Ng)),1,h_tsv);

%% Cleanup - Get rid of NaNs
iidf(isnan(iidf)) = 0;
lmax = length(iidf) - 1;
l = 0:lmax;

%% Determine wire pitch and layer assignment

Ach_wla = Ach_tier_gp; % Reduce the chip area by a factor of S when we're folding a design across S layers
[Ln pn pn_orig Nm] = wire_layer_assignment_alt(iidf,lmax,Ach_wla,chi,rho_m,epsr_d,Tclk,alpha_t);
[Cxc Ltot Cn] = calc_total_wiring_capacitance2(pn,Ln,Nm,iidf,epsr_d,gate_pitch);

%% Power estimates

eps0 = 8.854e-12; % (F/m)

Ilk = Ioff*(w_trans*1e6);
Cox = eps_ox*eps0*w_trans^2/tox;
Co = Cox; % Need to include parasitics for realistic estimate

Nt = N_trans_per_gate * Ng;
f = 1/Tclk;

Pdyn = 1/2*a*Co*Vdd^2*f*Nt;
Plk = (1-a)*Vdd*Ilk*Nt;
Pw = 1/2*a*Cxc*Vdd^2*f;

%% Repeater insertion

Ach_ri = Ach_tier_m2;
Ainv_min = gate_pitch^2*9; % assume 3:1 W/L for nmos, 3x that for pmos
rho_xcn = rho_m;

Co = N_trans_per_gate*Cox;

% construct some objects for the repeater insertion routine
wire.Ln = Ln;
wire.pn = pn;
wire.dielectric_epsr = epsr_d;
chip.gate_pitch = gate_pitch;

[iidf_rep h_vec k_vec Arep_used num_vec size_vec] = repeater_insertion_old_capfix(iidf,Ach_ri,Ainv_min,pn,Ln,Cn,rho_xcn,Ro,Co,gate_pitch,chip,wire);

Co_rep = Cox*size_vec;
Ilk_rep = Ilk*size_vec;
Plk_rep_vec = (1-a)*Vdd*Ilk_rep.*num_vec;
Plk_rep = sum(Plk_rep_vec);

Pdyn_rep_vec = 1/2*a*Vdd^2*f.*num_vec.*Co_rep;
Pdyn_rep = sum(Pdyn_rep_vec);

Prep = Pdyn_rep + Plk_rep;

Arep_used_mm2 = Arep_used*(1e3)^2;

%% Redo wiring now that we've changed Iidf
if(redo_wiring == 1)
    iidf_rewire = [0 iidf_rep]; % add zero-length value back in
    [Ln pn pn_orig Nm] = wire_layer_assignment_alt(iidf_rewire,lmax,Ach_wla,chi,rho_m,epsr_d,Tclk,alpha_t);
    [Cxc Ltot Cn] = calc_total_wiring_capacitance2(pn,Ln,Nm,iidf_rewire,epsr_d,gate_pitch);

    Pw = 1/2*a*Cxc*Vdd^2*f;
else
    iidf_rewire = 0;
end

