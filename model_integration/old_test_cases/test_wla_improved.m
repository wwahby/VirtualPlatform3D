clear all
close all
%% Chip descriptors

% Stack parameters
S = 1;

%% 65nm Merom, entire chip
% Ng = 291e6/4;
% Ach_mm2 = 143;
% gate_pitch = 700e-9; % average gate pitch (sqrt(A_core/Ngates))
% min_pitch = 220e-9; % actual contacted gate pitch
% fmax = 3.0e9;

%% 65nm Merom, single core, wrong area
% Ng = 10e6;
% Ach_mm2 = 16;
% gate_pitch = 630e-9; % average gate pitch (sqrt(A_core/Ngates))
% min_pitch = 220e-9; % actual contacted gate pitch
% fmax = 3.0e9;

%% 65nm Merom, single core
% Ng = 10e6;
% Ach_mm2 = 27;
% gate_pitch = 820e-9; % average gate pitch (sqrt(A_core/Ngates))
% min_pitch = 220e-9; % actual contacted gate pitch
% fmax = 3.0e9;

%% 32nm Sandy Bridge, entire chip
% Ng = 2.27e9/4;
% Ach_mm2 = 435;
% gate_pitch = 875e-9;
% min_pitch = 112.5e-9;
% fmax = 3.6e9;

%% 32nm Sandy Bridge, one core
% Ng = 107e6/4;
% Ach_mm2 = 20;
% gate_pitch = 435e-9;
% min_pitch = 112.5e-9;
% fmax = 3.6e9;

%% 22nm Ivy Bridge, entire chip
% Ng = 2.86e9/4;
% Ach_mm2 = 346.5;
% gate_pitch = 348e-9;
% min_pitch = 90e-9;
% fmax = 3.0e9;

%% 22nm Ivy Bridge, one core
Ng = 95e6/4;
Ach_mm2 = 11.5;
gate_pitch = 348e-9;
min_pitch = 90e-9;
fmax = 3.0e9;

%% Arbitrarily huge test case
% Ng = 1e9;
% Ach_mm2 = 100;
% gate_pitch = 100e-9;
% min_pitch = 100e-9;
% fmax = 3.0e9;

%% Carry on
Ach_m2 = Ach_mm2*1e-6;

% Tsv parameters
Atf_max = 0.00; % maximum allowable TSV area, as a fraction of total chip area
h_tsv_m = 10e-6;
AR_tsv = 20;

% Rent parameters
p = 0.6; % rent exponent
fo = 4; % avg fanout
alpha = fo/(fo+1); % input terminal fraction
%k = 3/alpha; %rent constant
k = 4;


%% Presize the chip and TSVs
Ns = Ng/S;
Lx = round(sqrt(Ns));

% Recalculate these to make sure everything is a nice integer
Ns = Lx^2;
Ng = Ns*S;

% Size the TSVs
h_tsv = ceil(h_tsv_m/gate_pitch);
w_tsv = ceil(h_tsv/AR_tsv);

%% Size the chip so we have a nicely divisible number of unit cells per side
%Nsp = floor( (1+Atf_max)*Ns );
Nsp = floor(Ns/(1-Atf_max));
Lxp = floor(sqrt(Nsp));
Nsp = Lxp^2;
Tp = ceil(w_tsv/sqrt(Atf_max));

slack = 0.2;
[Lxc Tc Nuc_1d gfrac_L gfrac_T] = find_LT_combination(Lxp,Tp,slack);

Nsc = Lxc^2;
Ngc = Nsc*S;

N_tsvs = Nuc_1d^2;
g_tsv = (Nuc_1d*w_tsv)^2; % number of gates displaced by TSVs
Atf_act = g_tsv/Nsc;

repstr1 = sprintf('Ng_nom: %.4g',Ng);
repstr2 = sprintf('Ng_act: %.4g',Ngc);
repstr3 = sprintf('Atf_act: %.4g',Atf_act);
disp(repstr1)
disp(repstr2)
disp(repstr3)

%% Calculate WLD
iidf_2d = calc_Iidf(alpha,k,p,round(sqrt(Ng)),1,h_tsv);
iidf_c = calc_Iidf_corrected(alpha,k,p,Lxc,S,h_tsv,Nuc_1d,w_tsv);
iidf_j = calc_Iidf(alpha,k,p,Lx,S,h_tsv);

%AA = calc_Nstart(Lx,S,r);
% Mt3dj = Mt_3d_joyner(Lx,S,r);
% Mt3dc = Mt_3d_corrected(Lx,S,r,Nuc_1d,w_tsv);
% Mt2dj = Mt_2d_joyner(Lx);
% Mt2dc = Mt2d_corrected(Lx, Nuc_1d, w_tsv);
% Nstart = calc_Nstart(Lx,S,r,g_tsv);
% Nnst = calc_Nnst(Lx,S,r,g_tsv);
% Nnsb = calc_Nnsb(Lx,S,r,g_tsv);
% Nc = calc_Nc(Mt3dc,Lx,S,r,g_tsv);
% h = calc_h(Lx, Nuc_1d, w_tsv);
% term4 = zeros(1,length(Mt2dj));
% term4(1:2*w_tsv+1) = N_tsvs*Mt_2d_joyner(w_tsv);
% iexp = calc_Iexp(alpha,k,p,Mt3dc,Lx,S,r,g_tsv);

%% Get rid of NaNs
iidf_2d(isnan(iidf_2d)) = 0;
lmax_2d = length(iidf_2d) - 1;
l2d = 0:lmax_2d;

iidf_c(isnan(iidf_c)) = 0;
lmax_c = length(iidf_c) - 1;
lc = 0:lmax_c;

iidf_j(isnan(iidf_j)) = 0;
lmax_j = length(iidf_j) - 1;
lj = 0:lmax_j;

%% Estimate total wiring capacitance

% fmax = 3.0e9;
chi = 2/3;
rho_m = 17.2e-9; % Cu
epsr_d = 2.9; % Low-k dielectric
Tclk = 1/fmax; % (s)
alpha_t = 1.1*6.2;

Ach = Ach_m2; % Force chip to use specified area 
Ach_gp = Ach_m2/gate_pitch^2;
Ach_mp = Ach_m2/min_pitch^2;
layers_per_tier = 1;
routing_efficiency_vec = [0.2 0.2 0.3 0.4];
routing_efficiency = 0.4;
routing_efficiency_vec = [0.4];
layer_area = Ach;


Beta = [0.25 0.9];
A_layer_max = layer_area*routing_efficiency*layers_per_tier;
Rc = 0;
%Rc = [1e2 1e1 1e0 1e-1];
rho_m_vec = [1 1 1] * rho_m;
%rho_m_vec = rho_m;

[Ln_vec pn_vec pn_orig_vec A_wires A_vias] = wla_improved_old(iidf_2d,gate_pitch,min_pitch,layers_per_tier,routing_efficiency,layer_area,rho_m_vec,epsr_d,alpha_t,Beta,Tclk,Rc);
[Ln2d pn2d pn2d_orig] = wire_layer_assignment_alt(iidf_2d,lmax_2d,Ach_gp,chi,rho_m,epsr_d,Tclk,alpha_t);
[Ln2dm pn2dm pn2d_origm] = wire_layer_assignment_alt(iidf_2d,lmax_2d,Ach_mp,chi,rho_m,epsr_d,Tclk,alpha_t);

Beta = 0.9;
Ro = 10e3;
Co = 1e-15;
%repeater_fraction = [0.5 0.3 0.2];
repeater_fraction = 0.5;
%% pack inputs
chip.iidf = iidf_2d;
chip.gate_pitch = gate_pitch;
chip.min_pitch = min_pitch;
wire.layers_per_tier = layers_per_tier;
wire.routing_efficiency = routing_efficiency_vec;
wire.layer_area = layer_area;
wire.resistivity = rho_m;
wire.dielectric_epsr = epsr_d;
wire.delay_constant = alpha_t;
wire.Beta = Beta;
chip.clock_period = Tclk;
wire.Rc = Rc;
chip.lengths = l2d;
chip.Ro = Ro;
chip.Co = Co;
wire.repeater_fraction = repeater_fraction;
chip.area_total = Ach_m2;

%% New routines

wire_norep = wla_improved(chip,wire);
wire_rep = wla_topdown_with_repeaters(chip,wire);

% [Ln_vec_td pn_vec_td A_wires_td A_vias_wiring_td A_vias_repeaters_td repeater_num_td repeater_size_td] = ...
%     wla_topdown_with_repeaters(...
%     iidf_2d,gate_pitch,min_pitch,layers_per_tier,routing_efficiency_vec,...
%     layer_area,rho_m,epsr_d,Beta,Tclk,Rc,Ro,Co,repeater_fraction);


% [Ln_vec pn_vec pn_orig_vec A_wires A_vias] = wla_improved(iidf_c,gate_pitch,min_pitch,layers_per_tier,routing_efficiency,layer_area/S,rho_m,epsr_d,alpha_t,Beta,Tclk,Rc);
% [Ln2d pn2d pn2d_orig] = wire_layer_assignment_alt(iidf_c,lmax_c,Ach_gp/S,chi,rho_m,epsr_d,Tclk,alpha_t);
% [Ln2dm pn2dm pn2d_origm] = wire_layer_assignment_alt(iidf_c,lmax_c,Ach_mp/S,chi,rho_m,epsr_d,Tclk,alpha_t);
%% Capacitance
[Cxc2d Ltot2d Cn2d] = calc_total_wiring_capacitance(pn2d,Ln2d,iidf_2d,l2d,epsr_d,gate_pitch);
[Cxc2dm Ltot2dm Cn2dm] = calc_total_wiring_capacitance(pn2dm,Ln2dm,iidf_2d,l2d,epsr_d,gate_pitch);
[Cxc Ltot Cno] = calc_total_wiring_capacitance(pn_vec,Ln_vec,iidf_2d,l2d,epsr_d,gate_pitch);
[Ctot Cnn] = calc_wiring_capacitance_from_area_old(pn_vec,layers_per_tier,A_wires,A_vias,epsr_d);

%% Size gates and see how far off we were
Vt = 0.2;
Vdd = 1.325;
epsr_ox = 3.9;
t_ox = 1.2e-9;
mobility_n = 1300; % (cm^2/Vs);
mobility_p = mobility_n/3;
mobility_av = 1/2*(mobility_n + mobility_p);

L_gate_min = 22e-9;
Wmin = L_gate_min;
logic_depth = 5;
Beta_g = 0.75;


Rnand = calc_Rnand(Vdd,Vt,epsr_ox,t_ox,mobility_av);
Cnand = calc_Cnand(L_gate_min,epsr_ox,t_ox);

l_av = get_average_wirelength(iidf_2d);
C_av = get_capacitance_from_length_old(l_av,Ln_vec,pn_vec,epsr_d,gate_pitch);
[W gate_pitch_new] = size_gates_nand2(Wmin,Rnand,Cnand,100*C_av,fo,logic_depth,Tclk,Beta_g);


%%
figure(1)
clf
plot(pn2d,'b')
hold on
plot(pn2dm,'g--')
plot(pn_vec/min_pitch,'r')
plot(wire_norep.pn/min_pitch,'m--')
plot(wire_rep.pn/min_pitch,'k--')
xlabel('Wiring tier')
ylabel('Wire pitch (MP)')
legend('Orig-GP','Orig-MP','Improved','newio-norep','newio-tdrep')
fixfigs(1,3,14,12)


figure(2)
clf
semilogy(A_wires./A_layer_max,'b');
hold on
semilogy(A_vias./A_layer_max,'r');
xlabel('Wiring tier')
ylabel('Area fraction')
legend('wires','vias')
fixfigs(2,3,14,12)

% figure(3)
% clf
% plot(Ln2d,'b')
% hold on
% plot(Ln2dm,'g--')
% plot(Ln_vec,'r')
% plot(wire_norep.Ln,'m--')
% plot(wire_rep.Ln,'k--')
% xlabel('Wiring tier')
% ylabel('Longest wire routed per tier')
% legend('Orig-GP','Orig-MP','Improved')
% fixfigs(3,3,14,12)

% figure(4)
% clf
% plot(Cn2d.*Ltot2d*gate_pitch/1e-9,'b')
% hold on
% plot(Cn2dm.*Ltot2dm*gate_pitch/1e-9,'g')
% plot(Cno.*Ltot.*gate_pitch/1e-9,'m')
% plot(Cnn/1e-9,'r')
% xlabel('Wiring tier')
% ylabel('Capacitance per tier (nF)')
% legend('Orig-GP','Orig-MP','Improved-oldC','Improved-newC')
% fixfigs(4,3,14,12)

figure(5)
clf
plot(pn2d*gate_pitch/1e-9,'b')
hold on
plot(pn2dm*min_pitch/1e-9,'g--')
plot(pn_vec/1e-9,'r')
plot(wire_norep.pn/1e-9,'m--')
plot(wire_rep.pn/1e-9,'k--')
xlabel('Wiring tier')
ylabel('Wire pitch (nm)')
legend('Orig-GP','Orig-MP','Improved','Location','nw')
fixfigs(5,3,14,12)

