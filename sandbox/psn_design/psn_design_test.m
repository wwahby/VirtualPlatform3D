% Add VP3D to path for simulation setup
addpath('../../model_integration/')

area = (sqrt(2)*1e-3)^2; %m2
core_power = 2e-3; % W

wire.rho_vec = 16.8e-9; % Cu
wire.permeability_rel = 1; % nonmagnetic
chip.area_per_layer_m2 = area;
chip.num_layers = 64;
chip.temperature = 300; % K
chip.Vdd = 1;

psn.decap_area_fraction = 0.1;
tsv.height_m = 50e-6;
tsv.barrier_thickness = 0e-9;
tsv.barrier_resistivity = 17.2e-9;

psn.noise_target = 0.1875;          % (V) Acceptable power supply noise
psn.noise_fraction = 0.10;
psn.decap_area_fraction = 0.1;      % (Ratio) - Fraction of chip area dedicated to decoupling capacitors
psn.power_tsv_width = 10e-6;        % (m) Width of TSVs used for power delivery. This may be larger than signal TSVs. Setting this to -1 will cause PSN module to use the signal TSV width.

% Power connection parameters
psn.Npads_1d = 100;                  % (-) Number of pads from one side of chip to the other
psn.Npads = psn.Npads_1d^2;         % (-) total number of power pads (does not include ground pads)
psn.Ngrid = 21*21;                  % (-) grid fineness
psn.pad_size = 1;                   % (-) Pad size is in terms of segment number
psn.segment_thickness = 1e-6;       % (m) Thickness of a power grid segment
psn.segment_width = 2e-6;           % (m) width of a power grid segment

% Package parameters
psn.package_resistance = 0.006;     % (Ohm) Resistance per pad on the package
psn.package_inductance = 0.5e-9;    % (H) Inductance per package pad

% Power TSV determination
psn.mismatch_tolerance = 0.01;      % (-) Allowable normalized deviation from noise target

power.density = core_power/area;

psn = determine_power_tsv_requirements(tsv,psn,power,wire,chip);