function [psn_max] = calc_psn(psn,power,chip,tsv,rho_m,mu_m,T)
% === INPUTS ===
% area      (m^2)   chip area
% Npad      (-)     Number of power pads -- does not include ground pads!
% rho_m     (Ohm*m) Resistivity of the power/gnd wires
% mu_m      (H/m)	Permeability of the power/gnd wires
% Nstrata   (-)     Number of chips in the 3D stack
% Vdd       (V)     Supply voltage
% Pdens     (W/m^2) Power density
% decap     (-)     Fraction of chip area devoted to decoupling capacitors
% h_tsv     (m)     TSV height
% w_tsv     (m)     TSV width
% pitch_tsv (m)     TSV pitch
% RPKG      (Ohm)   Package resistance - get this from literature or sims
% LPKG      (H)     Package inductance - get this from literature or sims
% Ngrid     (-)     Number of grid points (grid fineness)
% padsize   (-)     Power/gnd pad size in units of unit grid cells
% Tseg      (m)     Thickness of power network grid segment
% Wseg      (m)     Width of power network grid segment
% T         (deg C) Temperature
% ==
% === OUTPUTS ===
% psn       (W)     Power supply noise
% ==

%%% Some test inputs below
% %% Chip parameters
% area = 184e-6;        %chip area (m^2)
% Npad = 32*32;         %total number of power or ground pads
% rho_m = 1.68*10^-8;     %copper resistivity
% mu_m = 1.257e-6;      %copper permeability
% Nstrata = 2;
% 
% Vdd = 1.0; % (V)
% Pdens = 100; % (W/cm2)
% decap = 0.1; % (Ratio) - Fraction of chip area dedicated to decoupling capacitors
% 
% T = 25; % (deg C) Temperature
% 
% %% Unit Cell Constants
% Ngrid = 21*21;        %grid fineness
% padsize = 1;          % Pad size is in terms of segment number
% Tseg = 1e-6;          % Thickness of grid segment
% Wseg = 2e-6;          % width of grid segment
% 
% %% TSV constants
% h_tsv = 50e-6;         %TSV height
% w_tsv = 10e-6;         %TSV diameter
% pitch_tsv = 500e-6;
% 
% 
% %% Package parameters
% RPKG = 0.006; % (Ohm)
% LPKG = 0.5e-9; % (H)

%% Unpack inputs from objects
%psn_max = calc_psn(chip.area_per_layer_m2,Npad,rho_m,mu_m,S,Vdd,Pdens,decap,h_tsv_m,w_tsv_m,pitch_tsv,RPKG,LPKG,Ngrid,padsize,Tseg,Wseg,T);

area = chip.area_per_layer_m2;
Npad = psn.Npads;
Nstrata = chip.num_layers;
Vdd = chip.Vdd;
Pdens = power.density;
decap = psn.decap_area_fraction;
h_tsv = tsv.height_m;
w_tsv = tsv.width_m;
pitch_tsv = psn.pitch_tsv;
RPKG = psn.package_resistance;
LPKG = psn.package_inductance;
Ngrid = psn.Ngrid;
padsize = psn.pad_size;
Tseg = psn.segment_thickness;
Wseg = psn.segment_width;


%% Calculate useful quantities
Jch = Pdens/Vdd; % (A/cm^2) Current density
pp = sqrt(2*pitch_tsv^2); % [FIX] Not really sure what the point of this is

%% Determine Unit Cell parameters

[lcell,Ppitch,PGpitch,rpad,Rseg,Cd] = Ucell(area,Npad,Ngrid,padsize,rho_m,Tseg,Wseg,decap);
acell = lcell;

%% Determine Package Parameters
% [FIX] using arbitrary parameters for now
[RTSV,LTSV] = RL_TSV(h_tsv,w_tsv,rho_m,mu_m,pp,pitch_tsv);

%% Determine power supply noise

layer = Nstrata; % Just worry about top die (worst case scenario)
psn_max = psn_fast(Nstrata,layer,RTSV,LTSV,RPKG,LPKG,Cd,Jch,T,acell,Rseg,rpad);

% Print the inputs to psn_fast for debug purposes -- leave this commented
% out most of the time
% psnstr = sprintf('Nstrata_%d	layer_%d	RTSV_%d	LTSV_%d	RPKG_%d	LPKG_%d	Cd_%d	Jch_%d	T_%d	acell_%d	Rseg_%d	rpad_%d',Nstrata,layer,RTSV,LTSV,RPKG,LPKG,Cd,Jch,T,acell,Rseg,rpad);
% disp(psnstr)