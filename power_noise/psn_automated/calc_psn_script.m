% script to run Li's PSN model
% want this to be fully automated
%clear all
%% Chip parameters
area = (100e-3)^2;
Npad = 32*32;         %total number of power or ground pads
rho = 1.68*10^-8;     %copper resistivity
mu = 1.257e-6;      %copper permeability
Nstrata = 2;

Vdd = 1.0; % (V)
Pdens = 100; % (W/cm2)
decap = 0.1; % (Ratio) - Fraction of chip area dedicated to decoupling capacitors
Jch = Pdens/Vdd;
T = 25; % (deg C) Temperature

%% Unit Cell Constants
Ngrid = 21*21;        %grid fineness
padsize = 1;          % Pad size is in terms of segment number
Tseg = 1e-6;          % Thickness of grid segment
Wseg = 2e-6;          % width of grid segment

%% TSV constants
TH = 50e-6;         %TSV height
TD = 10e-6;         %TSV diameter
pitch = 500e-6;
pp = sqrt(2*pitch^2);

%% Package parameters
RPKG = 0.006; % (Ohm)
LPKG = 0.5e-9; % (H)

%% Determine Unit Cell parameters

[lcell,Ppitch,PGpitch,rpad,Rseg,Cd] = Ucell(area,Npad,Ngrid,padsize,rho,Tseg,Wseg,decap);
acell = lcell;

%% Determine Package Parameters
% [FIX] using arbitrary parameters for now
[RTSV,LTSV] = RL_TSV(TH,TD,rho,mu,pp,pitch);

%% Determine power supply noise

layer = Nstrata; % Just worry about top die (worst case scenario)
tic
psn1 = psn(Nstrata,layer,RTSV,LTSV,RPKG,LPKG,Cd,Jch,T,acell,Rseg,rpad)
time1 = toc

tic
psn2 = psn_fast(Nstrata,layer,RTSV,LTSV,RPKG,LPKG,Cd,Jch,T,acell,Rseg,rpad)
time2 = toc