%clear all
%close all
%% Model flags
use_joyner = 0; % 0=Use corrected distribution, 1=use original joyner distribution
redo_wiring = 0; % Redo wire layer assignment after repeater insertion? 0=no, 1=yes

%% Sweep parameters
% S_vec = [1 2 4];
% h_vec = linspace(10e-6,300e-6,2);

S_vec = [1 2 4 8];
%h_vec = linspace(10e-6,300e-6,11);
h_vec = logspace(-2,0,21)*300e-6;
h_vec = [1e-6 2e-6 h_vec];

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
Atf_max = 0.011; % maximum allowable TSV area, as a fraction of total chip area
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
fmax = 3e9;
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


%% Sweep tsv height and number of layers
nh = length(h_vec);
ns = length(S_vec);

Pw_mat = zeros(nh,ns);
Cxc_mat = zeros(nh,ns);

for hind = 1:nh
    for sind = 1:ns
        S = S_vec(sind);
        h_tsv_m = h_vec(hind);
        [ iidf l Ln pn pn_orig Cxc Ltot Cn Pdyn Plk Pw Prep Ng_act N_tsvs ] = gen_design(Ng,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);
        
        Pw_mat(hind,sind) = Pw;
        Cxc_mat(hind,sind) = Cxc;
    end
end

%%
hscale = 1e-6;
figure(1)
clf
set(gcf,'DefaultAxesColorOrder',[0 0 0; 1 0 0; 0 0 1; 0 1 0]); % color
%set(gcf,'DefaultAxesColorOrder',[0 0 0]); % b/w
set(gcf,'DefaultAxesLineStyleOrder','-|--|:|-.')
set(gcf,'DefaultLineLineWidth',2)
hold all
for sind = 1:ns
    plot(h_vec/hscale,Pw_mat(:,sind))
end

% logarithmic x axis
set(gca,'xscale','log'); xlim([1 300]);
ylim([25 55])

%legend('2D','2 tiers','4 tiers','8 tiers','location','se')
xlabel('Interlayer separation (microns)')
ylabel('Wiring power (W)')
fixfigs(1,3,12,12)


