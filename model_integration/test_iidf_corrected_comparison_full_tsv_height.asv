%clear all
%close all
%% Model flags
use_joyner = 0; % 0=Use corrected distribution, 1=use original joyner distribution
redo_wiring = 0; % Redo wire layer assignment after repeater insertion? 0=no, 1=yes

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
Atf_max = 0.11; % maximum allowable TSV area, as a fraction of total chip area
gate_pitch = 90e-9;
h_tsv_m_thin = 10e-6;
h_tsv_m_thick = 300e-6;
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




%% Corrected distribution
h_tsv_m = h_tsv_m_thin;

S=2; 
[ iidf_3d2c l_3d2c Ln_3d2c pn_3d2c pn_orig_3d2c Cxc_3d2c Ltot_3d2c Cn_3d2c Pdyn_3d2c Plk_3d2c Pw_3d2c Prep_3d2c Ng_act_3d2c N_tsvs_3d2c ] = gen_design(Ng,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);

S=4;
[ iidf_3d4c l_3d4c Ln_3d4c pn_3d4c pn_orig_3d4c Cxc_3d4c Ltot_3d4c Cn_3d4c Pdyn_3d4c Plk_3d4c Pw_3d4c Prep_3d4c Ng_act_3d4c N_tsvs_3d4c ] = gen_design(Ng,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);

S=1;
[ iidf_2dc l_2dc Ln_2dc pn_2dc pn_orig_2dc Cxc_2dc Ltot_2dc Cn_2dc Pdyn_2dc Plk_2dc Pw_2dc Prep_2dc Ng_act_2dc N_tsvs_2dc ] = gen_design(Ng_act_3d2c,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);

%% Joyner distribution
h_tsv_m = h_tsv_m_thick;

S=1;
[iidf_2dj l_2dj Ln_2dj pn_2dj pn_orig_2dj Cxc_2dj Ltot_2dj Cn_2dj Pdyn_2dj Plk_2dj Pw_2dj Prep_2dj Ng_act_2dj N_tsvs_2dj ] = gen_design(Ng_act_3d2c,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);

S=2;
[ iidf_3d2j l_3d2j Ln_3d2j pn_3d2j pn_orig_3d2j Cxc_3d2j Ltot_3d2j Cn_3d2j Pdyn_3d2j Plk_3d2j Pw_3d2j Prep_3d2j Ng_act_3d2j N_tsvs_3d2j ] = gen_design(Ng_act_3d2c,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);

S=4;
[ iidf_3d4j l_3d4j Ln_3d4j pn_3d4j pn_orig_3d4j Cxc_3d4j Ltot_3d4j Cn_3d4j Pdyn_3d4j Plk_3d4j Pw_3d4j Prep_3d4j Ng_act_3d4j N_tsvs_3d4j ] = gen_design(Ng_act_3d4c,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);

%% Calculate error
% ldiff1 = length(iidfc) - length(iidfj);
% ldiff2 = length(iidfj) - length(iidfc);
% iidfc = [iidfc zeros(1,ldiff2)];
% iidfj = [iidfj zeros(1,ldiff1)];
% 
% err_raw = abs(iidfc-iidfj);
% err_norm = err_raw./iidfj;

%% WLD comparison
lw = 2;

figure(1)
clf
loglog(iidf_2dj,'k','linewidth',lw);
hold on
%loglog(iidf_2dc,'k--','linewidth',lw,'linesmoothing','on');
loglog(iidf_3d2j,'b','linewidth',lw);
loglog(iidf_3d4j,'g','linewidth',lw);
loglog(iidf_3d2c,'r','linewidth',lw);
loglog(iidf_3d4c,'m','linewidth',lw);
xlim([1e4 1e5])
ylim([1e-1 1e2])

legend('2D','2 tier - 300um','4 tier - 300um','2 tier - 10um','4 tier - 10um')
ylabel('Number of interconnects')
xlabel('length (gate pitches)')
%title('Iidf')

%% Power comparison
figure(2)
clf
Pw2d = [Pw_2dj Pw_2dj];
Pw3d2 = [Pw_3d2j Pw_3d2c];
Pw3d4 = [Pw_3d4j Pw_3d4c];
    
h = [300 10];
plot(h,Pw2d,'k-o')
hold on
plot(h,Pw3d2,'k--s')
plot(h,Pw3d4,'k-.d')
legend('2D','2 tier','4 tier')
title('Wiring power')
%set(gca,'xscale','log')
%ylim([0 1.3*Pw_2dj])
xlim([10 300])
xlabel('TSV height')

%% Pitch comparison

figure(2)
clf
plot(pn_2dj,'k','linewidth',lw)
hold on
plot(pn_3d2j,'b','linewidth',lw)
plot(pn_3d4j,'g','linewidth',lw)
plot(pn_3d2c,'r','linewidth',lw)
plot(pn_3d4c,'m','linewidth',lw)

plot(pn_orig_2dj,'k--','linewidth',lw)
plot(pn_orig_3d2j,'b--','linewidth',lw)
plot(pn_orig_3d4j,'g--','linewidth',lw)
plot(pn_orig_3d2c,'r--','linewidth',lw)
plot(pn_orig_3d4c,'m--','linewidth',lw)

legend('2D','2 tier - 300um','4 tier - 300um','2 tier - 10um','4 tier - 10um')
xlabel('wiring tier')
ylabel('Wire pitch (gate pitches)')
%title('Wire pitch')

%% Longest wire routed
figure(3)
lw = 2;
clf
semilogy(Ln_2dj,'k','linewidth',lw)
hold on
semilogy(Ln_3d2j,'b','linewidth',lw)
semilogy(Ln_3d4j,'g','linewidth',lw)
semilogy(Ln_3d2c,'r','linewidth',lw)
semilogy(Ln_3d4c,'m','linewidth',lw)

legend('2D','2 tier - 300um','4 tier - 300um','2 tier - 10um','4 tier - 10um','location','se')

ylabel('longest wire routed (gate pitches)')
xlabel('Wiring tier')
title('Longest wire routed by metal layer')

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


%% Total capacitance
figure(5)
clf
% bar([Cxc_2dj Cxc_3d2j Cxc_3d4j Cxc_3d2c Cxc_3d4c])
% set(gca,'xticklabel',{'2D','2 tier - 300um','4 tier - 300um','2 tier - 10um','4 tier - 10um'})
C2d = [Cxc_2dj Cxc_2dj];
C3d2 = [Cxc_3d2j Cxc_3d2c];
C3d4 = [Cxc_3d4j Cxc_3d4c];
h = [300 10];
plot(h,C2d,'k-o')
hold on
plot(h,C3d2,'k--s')
plot(h,C3d4,'k-.d')
legend('2D','2 tier','4 tier')
title('Total wiring capacitance')
%set(gca,'xscale','log')
ylim([0 1.3*Cxc_2dj])
xlim([10 300])
xlabel('TSV height')

%% Another power figure


P2d = [Pw_2dj Prep_2dj];
P3d2 = [Pw_3d2c Prep_3d2c];
P3d4 = [Pw_3d4c Prep_3d4c];
gap = [0 0];

Pplot = [P2d ; gap ;  P3d2 ; gap ; P3d4];
Pplot = [P2d;  P3d2 ; P3d4];
Pplot = [P2d ; P3d4];

figure(6)
clf
barh = bar(Pplot,0.8,'stack');
set(barh,{'FaceColor'},{'b';'r'})
%set(gca,'xtick',[1 3 5])
set(gca,'xticklabel',{'1','4','4'})
legend('Wiring','Repeaters','location','ne')
xlim([0.5 2.5])
xlabel('# tiers')
ylabel('Power (W)')
fixfigs(6,3,14,12)
%%
fixfigs(1:5,2.5,12,12)

