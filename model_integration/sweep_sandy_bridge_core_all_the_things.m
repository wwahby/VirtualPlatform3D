close all
clear all

%% Simulation parameters
simulation.skip_psn_loops = 1; % Skip PSN TSV homing for faster debug
simulation.run_transient_psn = 0; % Run transient PSN in addition to static PSN
simulation.skip_thermal = 1; % Skip thermal analysis for faster debug
simulation.iterate_temperature = 0; % rerun WLA/Leakage until temperature converges
simulation.ignore_leakage = 0;

simulation.use_joyner = 0;
simulation.redo_wiring_after_repeaters = 0;
simulation.topdown_WLARI = 1; % Use topdown simultaneous WLA and RI (0 = use standard bottom-up optimal WLA, followed by one pass of RI)
simulation.separate_wiring_tiers = 1; % 1 = Each logic plane will have its own wiring tiers between it and the next logic plane
                                      % 0 = All metal layers for entire device will be routed on top of entire 3D stack
                                      
simulation.draw_thermal_map = 0; % Plot thermal profile of each chip
simulation.print_thermal_data = 0; % Output max temp in each layer to console

simulation.wla_max_attempts = 30; % 15 is default
simulation.wla_min_bot_fill_factor = 0.97; % 0.97 is default
simulation.wla_min_top_fill_factor = 0.01; % 0.01 is default

simulation.freq_binsearch = 0;
simulation.freq_binsearch_initial_guess = 1e9;
simulation.freq_binsearch_target = 90;
simulation.freq_binsearch_raw_tol = 0.25;
simulation.freq_binsearch_max_gens = 12;
simulation.freq_ceiling = 0;%3.5e9;

simulation.power_binsearch = 0;
simulation.power_binsearch_target = 20;
simulation.power_binsearch_raw_tol = 0.01;

simulation.heat_transfer_binsearch = 0;
simulation.heat_transfer_binsearch_temp_target = 90;
simulation.heat_transfer_binsearch_temp_raw_tol = 0.25;
simulation.heat_transfer_binsearch_max_gens = 10;

simulation.force_power = 0;
simulation.insanity_temperature = 1e3; % (C) 

simulation.ignore_repeater_area = 0; % 0 = restrict area, 1 = ignore area limit

%% Typical Rent Exponents
rent_exp_logic = 0.6;
rent_exp_mem = 0.4;
rent_exp_gpu = 0.55;

%% Logic core parameters
design.compression_factor = 1; % linear scaling factor. 1 = actual 32nm design, 4.57 = equivalent 7nm SB
design.Ng_core = 86e6/4; %86M transistors, assume 2in NAND -> /4 to get total NAND gates
design.Ach_mm2 = 18.5/(design.compression_factor^2);
design.gate_pitch = 465e-9*2/design.compression_factor;
design.min_pitch = 112.5e-9/design.compression_factor;
design.fmax = 3.5e9;
design.w_trans = 32e-9/design.compression_factor;
design.rent_exp = rent_exp_logic;

%% Thermal parameters
% Thermal resistances of known heat sinks (From Yue's paper)
r_air = 1/1.825; %K/W for a 1cm^2 HS % alt, 0.6
r_water = 1/4.63; %K/W for a 1cm^2 HS
A_hs = (1e-2)^2; % 1 cm^2

% Heat transfer coefficients
h_air = 1/(r_air*A_hs);
h_water = 1/(r_water*A_hs);

%% Metal resistivities
rho_ag = 15.9e-9;
rho_cu = 17.2e-9;
rho_au = 24.4e-9;
rho_al = 26.5e-9;
rho_w = 56.0e-9;
rho_ni = 69.9e-9;
rho_pt = 106e-9;
%rho_co_al = 0e-9;
rho_all_mets = [rho_ag rho_cu rho_au rho_al rho_w rho_ni rho_pt];

%% Thermal Conductivities
k_tim = 3;       % TIM
k_chip = 149;      % typically silicon
k_underfill = 0.9; % underfill
k_wires = 400;     % copper wires
k_copper = 400;    % Copper
k_tungsten = 173;  % Tungsten:
k_ild = 1.38;      % ILD
k_microbumps = 60; % Microbumps
k_interposer = 149; % interposer 
k_air = 0.024;      % Air

%% Sweep settings
tiers = [1 2 3 4];
num_gates_vec = design.Ng_core;
thicknesses = [0.1 100]*1e-6;
tsv_aspect_ratio = 10;
force_thickness = 1;
rel_permittivities = [3];
frequencies = design.fmax;
heat_fluxes = [ h_air];% h_water h_water];
temperature_targets = [70]; % temperature target to be used for each heat flux condition
cooling_configs = {'up'};%, 'down', 'down_all'}; % Location of heat sink in each heat flux condition
thermal_conductivities = 0.3;
decap_ratios = [0.1]; % Fraction of die area used for decoupling capacitors
wire_resistivities = [rho_cu, rho_w]; %[10:10:60]*1e-9;
wire_thermal_conductivities = [k_copper, k_tungsten];
wire_material_flags = {'00'}; % binary strings. bit1 = use_graphene, bit0 = use alt_em_mat
scaling_factors = [32/32]; % 32/22 32/14 32/10 32/7 32/5];
node_labels = {'32nm'}; %{'32nm', '22nm', '14nm', '10nm', '7nm', '5nm'}; % labels for plots involving scaling factors
% scaling_factors = 1; 
% node_labels = {'32nm'}; % labels for plots involving scaling factors
barrier_thicknesses = [0e-9];
barrier_resistivities = [1000e-9];
power_forced_vec = linspace(0,25,21); % (W) If simulation.force_power is 1, this forces power consumption to these values during sweep for thermal purposes
power_tsv_width = -1; % set this to -1 to use the same dimensions as signal TSVs (signal TSV dimensions are determined automatically in the interconnect module, and will typically be on the smallish end of things)

%This doesn't get swept, but rather, if you're doing a scaling sweep, you can input extra Vdd values here to have each scaled node use a different Vdd. If Vdd is constant (or stops scaling after a certain node) you can just have a single entry (or only the first few entries until it stops changing)
Vdd_vec = [ 1.25, 1.0, 0.95, 0.90, 0.85, 0.80]; % Vdd used at each **scaling node**.


%% Run the parameter sweep
%sweep_data = sweep_design(design, sweep, simulation);
sweep_design

%% plot stuff!
% plot_sweep_data( sweep, sweep_data, simulation )
% plot_sweep_data % just run as script
% plot.plot_heat_transfer_coeff_vs_tiers
% plot.plot_temp_vs_forced_power
% plot.plot_epc_vs_tiers
% plot.plot_freq_vs_scaling
% plot.plot_freq_vs_tiers_and_cooling
% plot.plot_power_tsvs_vs_scaling
% plot.plot_temp_underfill_vs_tiers
% plot.plot_metal_levs_and_power_vs_permittivity
% plot.plot_everything_vs_resistivity_and_tiers
% plot.plot_psn_test_stuff
% plot.plot_xc_power_vs_num_gates
plot.tsv_vs_m3d_power