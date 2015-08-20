close all
clear all

%% Simulation parameters
simulation.skip_psn_loops = 1; % Skip PSN TSV homing for faster debug
simulation.skip_thermal = 0; % Skip thermal analysis for faster debug

simulation.use_joyner = 0;
simulation.redo_wiring_after_repeaters = 0;
simulation.topdown_WLARI = 1; % Use topdown simultaneous WLA and RI (0 = use standard bottom-up optimal WLA, followed by one pass of RI)
simulation.separate_wiring_tiers = 1; % 1 = Each logic plane will have its own wiring tiers between it and the next logic plane
                                      % 0 = All metal layers for entire device will be routed on top of entire 3D stack
                                      
simulation.draw_thermal_map = 0; % Plot thermal profile of each chip
simulation.print_thermal_data = 0; % Output max temp in each layer to console

simulation.wla_max_attempts = 15; % 15 is default
simulation.wla_min_bot_fill_factor = 0.91; % 0.97 is default
simulation.wla_min_top_fill_factor = 0.01; % 0.01 is default

simulation.freq_binsearch = 0;
simulation.freq_binsearch_initial_guess = 1e10;
simulation.freq_binsearch_target = 90;
simulation.freq_binsearch_raw_tol = 0.25;
simulation.freq_binsearch_max_gens = 10;

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
design.Vdd = 1.25;
design.rent_exp = rent_exp_logic;

%% Thermal parameters
r_air = 1/1.825; %K/W for a 1cm^2 HS % alt, 0.6
r_water = 1/4.63; %K/W for a 1cm^2 HS
A_hs = (1e-2)^2; % 1 cm^2

h_air = 1/(r_air*A_hs);
h_water = 1/(r_water*A_hs);

%% Metal resistivities
rho_ag = 15.9e-9;
rho_au = 24.4e-9;
rho_cu = 17.2e-9;
rho_w = 56.0e-9;
rho_ni = 69.9e-9;
%rho_co_al = 0e-9;
rho_al = 26.5e-9;
rho_all_mets = [rho_ag rho_cu rho_au rho_al rho_w rho_ni];


%% Sweep settings
sweep.tiers = [1 2 3 4];
sweep.thicknesses = 10e-6;
sweep.force_thickness = 1;
sweep.rel_permittivities = [3.0];
sweep.frequencies = [1 2 4 6 8 10]*1e9;%fmax_core; % if simulation.freq_binsearch is set, (1) is min freq and (2) is max freq
sweep.heat_fluxes = [ h_air];
sweep.decap_ratios = [0.1];%[0.01 0.1 1];
sweep.wire_resistivities = [rho_cu];
sweep.wire_material_flags = {'00'}; % binary strings. bit1 = use_graphene, bit0 = use alt_em_mat
sweep.scaling_factors = [1];

%% Run the parameter sweep
t_sweep_start = cputime;
sweep_data = sweep_design(design, sweep, simulation);
t_sweep_stop = cputime;
fprintf('\nTotal time elapsed for parameter sweep: %.3g seconds\t(%.3g minutes)\n\n',(t_sweep_stop-t_sweep_start),(t_sweep_stop-t_sweep_start)/60);

%% plot stuff!
plot_sweep_data( sweep, sweep_data, simulation )


