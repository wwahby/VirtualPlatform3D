%clear all
%close all
%% Chip descriptors

% Stack parameters
Ng = 1.0e9;
Ach_mm2 = 100;
Ach_m2 = Ach_mm2*1e-6;

% gate parameters
eps_ox = 25; % HfO2
tox = 1e-9;
w_trans = 32e-9;
N_trans_per_gate = 4;
Ioff = 10e-9; %(A/um)

% Tsv parameters
Atf_max = 0.01; % maximum allowable TSV area, as a fraction of total chip area
gate_pitch = 90e-9;
h_tsv_m = 300e-6;
AR_tsv = 20;

% Rent parameters
p = 0.6; % rent exponent
fo = 4; % avg fanout
alpha = fo/(fo+1); % input terminal fraction
k = 3/alpha; %rent constant

% Wiring parameters
fmax = 1e9;
chi = 2/3;
rho_m = 17.2e-9; % Cu
epsr_d = 3.0; % Low-k dielectric
Tclk = 1/fmax; % (s)
alpha_t = 1.1*6.2;

% Repeater parameters
Ro = 1e3; % (Ohm) Gate output resistance

% Power parameters
a = 0.1; % logic activity factor
Vdd = 1.0; % (V)

% Model parameters
redo_wiring = 0;


%% Corrected distribution
use_joyner = 0;

S=2; 
[ iidf_3d2c l_3d2c Ln_3d2c pn_3d2c pn_orig_3d2c Cxc_3d2c Ltot_3d2c Cn_3d2c Pdyn_3d2c Plk_3d2c Pw_3d2c Prep_3d2c Ng_act_3d2c N_tsvs_3d2c ] = gen_design(Ng,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);

S=4;
[ iidf_3d4c l_3d4c Ln_3d4c pn_3d4c pn_orig_3d4c Cxc_3d4c Ltot_3d4c Cn_3d4c Pdyn_3d4c Plk_3d4c Pw_3d4c Prep_3d4c Ng_act_3d4c N_tsvs_3d4c ] = gen_design(Ng,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);

S=1;
[ iidf_2dc l_2dc Ln_2dc pn_2dc pn_orig_2dc Cxc_2dc Ltot_2dc Cn_2dc Pdyn_2dc Plk_2dc Pw_2dc Prep_2dc Ng_act_2dc N_tsvs_2dc ] = gen_design(Ng_act_3d2c,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);

%% Joyner distribution
use_joyner = 1;

S=1;
[iidf_2dj l_2dj Ln_2dj pn_2dj pn_orig_2dj Cxc_2dj Ltot_2dj Cn_2dj Pdyn_2dj Plk_2dj Pw_2dj Prep_2dj Ng_act_2dj N_tsvs_2dj ] = gen_design(Ng_act_3d2c,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);

S=2;
[ iidf_3d2j l_3d2j Ln_3d2j pn_3d2j pn_orig_3d2j Cxc_3d2j Ltot_3d2j Cn_3d2j Pdyn_3d2j Plk_3d2j Pw_3d2j Prep_3d2j Ng_act_3d2j N_tsvs_3d2j ] = gen_design(Ng_act_3d2c,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);

S=4;
[ iidf_3d4j l_3d4j Ln_3d4j pn_3d4j pn_orig_3d4j Cxc_3d4j Ltot_3d4j Cn_3d4j Pdyn_3d4j Plk_3d4j Pw_3d4j Prep_3d4j Ng_act_3d4j N_tsvs_3d4j ] = gen_design(Ng_act_3d4c,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);

%% Calculate error
[err_raw_2d err_norm_2d] = calc_err(iidf_2dj,iidf_2dc);
[err_raw_3d2 err_norm_3d2] = calc_err(iidf_3d2j,iidf_3d2c);
[err_raw_3d4 err_norm_3d4] = calc_err(iidf_3d4j,iidf_3d4c);


% ldiff1 = length(iidfc) - length(iidfj);
% ldiff2 = length(iidfj) - length(iidfc);
% iidfc = [iidfc zeros(1,ldiff2)];
% iidfj = [iidfj zeros(1,ldiff1)];
% 
% err_raw = abs(iidfc-iidfj);
% err_norm = err_raw./iidfj;

%% WLD comparison
figure(1)
clf
loglog(iidf_2dj,'k');
hold on
loglog(iidf_2dc,'k--');
loglog(iidf_3d2j,'b');
loglog(iidf_3d4j,'g');
loglog(iidf_3d2c,'r');
loglog(iidf_3d4c,'m');
title('Iidf')
fixfigs(1,2,12,12)

%% Power comparison
figure(2)
clf
P = [Pw_2dj Pw_3d2j Pw_3d4j Pw_3d2c Pw_3d4c];
bar(P)
title('Wiring Power')

%% Pitch comparison
lw = 2;
figure(2)
clf
plot(pn_2dj,'k')
hold on
plot(pn_3d2j,'b')
plot(pn_3d4j,'g')
plot(pn_3d2c,'r')
plot(pn_3d4c,'m')

plot(pn_orig_2dj,'k--')
plot(pn_orig_3d2j,'b--')
plot(pn_orig_3d4j,'g--')
plot(pn_orig_3d2c,'r--')
plot(pn_orig_3d4c,'m--')

legend('2D','2 tier - J','4 tier - J','2 tier - new','4 tier - new')
title('Wire pitch')
fixfigs(2,2,12,12)

%% Longest wire routed
figure(3)
lw = 2;
clf
semilogy(Ln_2dj,'k')
hold on
semilogy(Ln_3d2j,'b')
semilogy(Ln_3d4j,'g')
semilogy(Ln_3d2c,'r')
semilogy(Ln_3d4c,'m')

legend('2D','2 tier - J','4 tier - J','2 tier - new','4 tier - new','location','se')
title('longest wire routed per layer')
fixfigs(3,2,12,12)

%% Capacitance per layer
figure(4)
clf
plot(Cn_2dj,'k')
hold on
plot(Cn_3d2j,'b')
plot(Cn_3d4j,'g')
plot(Cn_3d2c,'r')
plot(Cn_3d4c,'m')
title('Capacitance per layer')
fixfigs(4,2,12,12)

%% Total capacitance
figure(5)
clf
bar([Cxc_2dj Cxc_3d2j Cxc_3d4j Cxc_3d2c Cxc_3d4c])
title('Total wiring capacitance')

%% Raw error
figure(6)
clf
loglog(err_raw_2d)
hold on
loglog(err_raw_3d2,'r')
loglog(err_raw_3d4,'g')

%% Normalized error
figure(7)
%clf
%loglog(err_norm_2d)
%loglog(100*err_norm_3d4,'b')
%hold on
%semilogx(100*err_norm_3d2,'r')
xlabel('length (gate pitches)')
ylabel('Deviation from original distribution (%)')
fixfigs(7,2,12,12)
ylim([0 20])



