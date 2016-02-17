function [max_noise, max_noise_time, time_mat, voltage_mat] = run_transient_power_analysis(Vdd, die_thickness, chip_width, chip_height, tsv, resistivity_bulk, power_per_die, decap_frac_per_die, temperature_K, heat )

%the testbench for PDN; set the systimatical parameters
%this file lists all the parameters needed
%spparms('spumoni', 0);
chip.die = die_thickness; %die thickness
chip.bump = 25e-6; %bump thickness between two dice
if (die_thickness < heat.monolithic_max_thickness)
    chip.bump = heat.monolithic_intertier_bond_thickness;
else
    chip.bump = heat.bump_thickness;
end

%%%%%%%%%%%%%%%%%%for TSV domain%%%%%%%%%%%%%%%%%%%%%
TSV.d = tsv.width_m; %TSV diameter
TSV.px = tsv.pitch_m; %TSV X pitch
TSV.py = tsv.pitch_m; %TSV Y pitch

% TSV.map = [3.0e-3 3.0e-3 1.4e-3 5.6e-3 200e-6 200e-6
%            7.3e-3 3.0e-3 1.4e-3 5.6e-3 200e-6 200e-6];
TSV.map = [];
%tsv subblock pitch needs to be 1/2; 1/4; 1/8 ..... of the global pitch
%this is the block which uses different TSV pitch
%the non-uniform meshing is based on the map
%the TSV.map is not yet implemented; SKIP this parameter
TSV.model_level = 1; %TSV model level, 2^TSV_model meshes between P and G pad
       
TSV.P = 0; %Power pad number counter
TSV.G = 0; %Ground pad number counter will be used in mesh function

electron_mfp = 39e-9; % (m) Mean free path of electrons in copper
specularity_coeff = 0.55;
reflection_coeff = 0.43;
barrier_thickness = 0;
rho_barrier = resistivity_bulk;
[R_tsv, rho_tsv, R_tsv_met, R_tsv_barrier] = xcm.calc_cu_wire_resistance_size_dependent(resistivity_bulk,TSV.d,TSV.d,barrier_thickness,rho_barrier,die_thickness,electron_mfp,specularity_coeff,reflection_coeff,temperature_K);

%[R_tsv,L_tsv] = RL_TSV(H,D,barrier_thickness,rho_barrier,rho,mu,TSVpitch,Ppitch,temperature_K);

alpha = 0.00393;
TSV.rou = rho_tsv; %17.1e-9*(1+alpha*(100-20)); % 17.1e-9 is the resistivity of copper @20 centigrade
TSV.R = R_tsv;
TSV.L = 0.25*1.256e-6/pi*(chip.die+chip.bump)*log(TSV.px/TSV.d);
%%%%%%%%%%%%%%%%%for Chip parameters%%%%%%%%%%%%%%%%%
chip.Xsize = chip_width; %chip x dimension; y dimension below
chip.Ysize = chip_height;

chip.wire_thick = 5e-6; %wire thickness
%chip.wire_p = 17.1e-9*(1+alpha*(100-20)); % wire resistivity
chip.wire_ar = 1.5; %aspect ratio = thickness/width (wire)
chip.wire_width = chip.wire_thick/chip.wire_ar;

[R_wire, rho_wire, R_wire_met, R_wire_barrier] = xcm.calc_cu_wire_resistance_size_dependent(resistivity_bulk,chip.wire_width,chip.wire_thick,barrier_thickness,rho_barrier,die_thickness,electron_mfp,specularity_coeff,reflection_coeff,temperature_K);
chip.wire_p = rho_wire;
chip.Rs = chip.wire_p/chip.wire_thick; %sheet resistance

chip.N = length(power_per_die); %die number
chip.cap_per = decap_frac_per_die;
chip.cap_th = 0.9e-9; %capacitance effective thicknee (used for capacitance value calculation)
chip.c_gate = 3.9*8.85e-12/chip.cap_th;

chip.power = power_per_die;% 74.49]; %74.49 74.49];  %total power dissipation of each die
%chip.power = [74.49 74.49]; %74.49 74.49];  %total power dissipation of each die

chip.map = zeros(length(power_per_die), 6);
for tier = 1:length(power_per_die)
    chip.map(tier, :) = [0 0 chip_width chip_height power_per_die(tier) decap_frac_per_die(tier)];
end

%blk_num is for splitting the power maps of each die
% format: xl corner, yb corner, x size, y size, power, decap %
chip.map(:,6) = chip.map(:,6)*chip.c_gate;
%chip.map(13:24,5) = chip.map(13:24,5);

chip.blk_num = ones(1,chip.N);

%%%%%%%%%%%%%%%%%for system level parameters%%%%%%%%%%%%
system.tran = 0; % transient analysis or not;
system.map = 0; % whether to draw the gifs along simulations
system.range = [0.5 1.3]; %the range for plotting gifs.
system.draw = 1; %whether to draw the power map for transient analysis
system.display = 1; %display the maximum noise of the chip
system.drawP = 1; %display the current requirement map
system.drawC = 0; %display the decap distribution
system.L = 0.1e-9; %package inductance for each pad
system.R = 0.0107; %package resistance for each pad
system.Vdd = Vdd; %VDD value for each pad

system.T = 50e-9; %transient simulation time 
system.Tr = 1e-9; %rise time
system.dt = 0.1e-9; %simulation step, TR scheme is stable for 0.1e-9, which is pretty optimized
system.ratio = 5; %the ratio between the maximum time step and minimum one
system.start_ratio = 1; %the start power of processor die; can optimize to generate more complicated power excitation
system.step = (1 - system.start_ratio)/(system.Tr/system.dt); %for each step, the increasement of the power excitation for processor die

[max_noise, max_noise_time, time_mat, voltage_mat] = pdnt.power_noise_sim(chip, system, TSV);
%main simulation function
