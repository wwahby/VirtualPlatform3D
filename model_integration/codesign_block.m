function [chip power tsv wire repeater psn] = codesign_block(chip,tsv,gate,transistor,wire,heat,psn,simulation)
time_start = cputime;
% Power and signal codesign
%% constants
eps0 = 8.854e-12; % (F/m) vacuum permittivity
mu0 = 4*pi*1e-7; % (H/m) Vacuum permeability

%% TSV number determination
% Inputs:
%   Rent parameters
%   Number of logic gates
%   Number of layers
disp(' ')
disp('Estimating TSV requirements...')
[nt_max nt_tot nt_to nt_through Tacmat] = xcm.estimate_tsvs_required(chip.num_gates,chip.num_layers,chip.rent_k,chip.rent_p,chip.alpha);

%% TSV Sizing
% Inputs:
%   Area available
%   Max area for TSVs
%   TSV aspect ratio
disp('Sizing TSVs...')
[w_tsv_m h_tsv_m] = xcm.size_tsvs(chip.area_total/chip.num_layers, tsv.max_area_fraction, nt_max, tsv.aspect_ratio );
%h_tsv_gp = round(h_tsv_gp);

if (simulation.force_thickness == 1)
    h_tsv_m = chip.thickness_nominal;
end
h_tsv_gp = ceil(h_tsv_m/chip.gate_pitch);
w_tsv_gp = ceil(w_tsv_m/chip.gate_pitch);

tsv.width_m = w_tsv_m;
tsv.height_m = h_tsv_m;
tsv.width_gp = w_tsv_gp;
tsv.height_gp = h_tsv_gp;
tsv.per_layer = nt_tot;
tsv.max_tsvs_per_layer = nt_max;

if(tsv.height_m > 0)
    chip.thickness = tsv.height_m;
else
    chip.thickness = chip.thickness_nominal;
end
%% System determination
% Run WLD + WLA + RI to get power estimate
disp('Generating system...')
%[ iidf l Ln pn pn_orig Cxc Ltot Cn Pdyn Plk Pw Prep Ng_act N_tsvs T_tsvs Atf_act ] = gen_design_old(Ng,alpha,k,p,S,h_tsv_m,Atf_max,AR_tsv,Ach_m2,chi,rho_m,epsr_d,Tclk,alpha_t,gate_pitch,w_trans,eps_ox,tox,N_trans_per_gate,a,Ioff,Vdd,Ro,use_joyner,redo_wiring);
[chip power wire repeater tsv] = xcm.gen_design(chip,tsv,gate,transistor,wire,simulation);

%Pdens = power.total/chip.area_per_layer_m2;
power.density = power.total/chip.area_total;

%% Thermal module -- Find actual system temperature

%%%%%%%%%%%%%%%%%%%%%%%geometry information of the chip%%%%%%%%%%%%%%%%%%%
    die.N = chip.num_layers;
    die.model = 3; % for each die, how many layers we model
    
    % flip chip package; order:
    %heatsink->TIM->CHIP_BULK1->METAL->BONDING->CHIP_BULK2-> ...
    %               CHIP_BULKN->METAL->MICRO-BUMPS->INTERPOSER    
    thick.bump = 40e-6;  %micro-bump thickness; between second die and interposer
    thick.tim = 5e-6; %tim thickness; between the chip and heatsink
    thick.under = 5e-6; %underfill bonding thickness; between two dies
    thick.inter = 200e-6; %interposer thickness
    thick.die = chip.thickness; %die thickness
    thick.ild = sum(wire.pn); %metal layer thickness
    
    layer_area = chip.area_per_layer_m2;
    side_length = sqrt(layer_area);
    grid_factor = 50;
    chip_therm.Xsize = side_length; %x dimension of chip
    chip_therm.Ysize = side_length; %y dimension of chip
    chip_therm.Xgrid = chip_therm.Xsize/grid_factor; %x grid size of chip
    chip_therm.Ygrid = chip_therm.Ysize/grid_factor; %y grid size of chip
    
    %for the interposer dimension
    pack.Xsize = 3.5*side_length;
    pack.Ysize = 3.5*side_length;
    pack.Xgrid = pack.Xsize/grid_factor;
    pack.Ygrid = pack.Ysize/grid_factor;
    
    %assumed TSV starting from top metal layer of a top die
    %            to the first metal layer of a bottom die
    %Thus the TSV passes through bonding layer & bulk of a botom die
    %TSV geometry 
    tsv_therm.d = tsv.width_m; % tsv diameter including the liner thickness
    tsv_therm.liner = tsv.width_m/10; %liner thickness
    tsv_therm.px = tsv.pitch_m; %x direction pitch
    tsv_therm.py = tsv.pitch_m; %y direction pitch
    tsv_therm.Nx = round(sqrt(tsv.num)); %x direction number
    tsv_therm.Ny = round(sqrt(tsv.num)); %y direction number

    %Bump geometry
    bump.d = 20e-6;
    bump.px = 100e-6;
    bump.py = 100e-6;
    bump.Nx = 40;
    bump.Ny = 40;
    
    portion = 0.5; %the metal portion in the ILD layers
    %This is used for equivalent thermal resistance calculation of metal
    %layer
    
    power_therm_vec = ones(1,chip.num_layers)*power.total/chip.num_layers;  %power dissipation of each die
    % from top to bottom; unit: watt
    
    granularity = 20;
    % thermal map, number of color used
    
    draw = simulation.draw_thermal_map;
    drawP = simulation.draw_thermal_map;
    % whether to draw the thermal map; 1 yes; 0 no
    
    displayT = simulation.print_thermal_data;
    % print the temperature information
%%%%%%%%%%%%%%%%%%%%%%%%finish geometry information%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%boundary condition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %the heat transfer coefficient
    % r = 1/(hA); A is the size of top surface area
    % the cooling capability of the top heatsink; 20000, 1cm*1cm, means:
    % 0.5 W/K
%     h.up = 20000;
%     
%     h.down = 5; % the cooling of bottom surface; 
%     %(only the area with the same size of chip;)
%     %microfluidic is assumed to be as large as chip in the interposer
%     
%     h.side = 5;
%     % side surface cooling, usually near adiabatic
%     
%     h.d = 5;
%     %the cooling of the bottom surface except for the MFHS area
%     
%     h.Ta = 298;
    %the ambient temperature
%%%%%%%%%%%%%%%%%%%%%%%%%finish boundary condition%%%%%%%%%%%%%%%%%%%%%%%%%

    chip_side = sqrt(chip.area_total/chip.num_layers);
    power_per_layer = power.total/chip.num_layers;
    

    %power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....
    map_row = [0     0     chip_side    chip_side     power_per_layer];
    map = zeros(chip.num_layers,5);
    blk_num = zeros(chip.num_layers,1);
    for i =1:chip.num_layers
        map(i,:) = map_row;
        blk_num(i,1) = 1;
    end
    %blk_num is for splitting the power maps of each die

    chip.temperature_vec = thermal.ThermSim( die, thick, chip_therm, pack, ...
              tsv_therm, bump, portion, power_therm_vec, ...
              map, blk_num, granularity, draw, drawP, heat, displayT);
          
    chip.temperature = max(chip.temperature_vec);


%% ============== END THERMAL MODULE ================

%% Power noise estimation
% Inputs:
%   Total power (from previous step)
%   TSV geometry from TSV Sizing step
disp('Evaluating power supply noise...')
psn.pitch_tsv = tsv.pitch_m*10; % [FIX] Power TSVs aren't going to be on the same pitch as signal TSVs
psn_mismatch_tolerance = psn.mismatch_tolerance; % 5% mismatch is OK
psn_iterations = 1;

rho_m = wire.resistivity(1); % use top layer
mu_m = wire.permeability_rel*mu0;

psn_max = power_noise.calc_psn(psn,power,chip,tsv,rho_m,mu_m,chip.temperature);
psn.noise = psn_max;
mismatch_norm = psn.noise/psn.noise_target;

dispstr = sprintf('\tpsn_runs: %d\tNpads: %d\tpsn_target: %d\tpsn_max: %d\tmismatch_norm: %.3g',psn_iterations,psn.Npads, psn.noise_target, psn_max,mismatch_norm);	
disp(dispstr)


%% Home in on the best number of power TSVs to use to meet the target
psn_iterations_max = 20;
while ( (abs(mismatch_norm-1) > psn_mismatch_tolerance) && (psn_iterations < psn_iterations_max) && (simulation.skip_psn_loops == 0))
    % more pads -> less noise, and vice versa
    psn.Npads_1d = round(sqrt(psn.Npads*mismatch_norm));
    psn.Npads = psn.Npads_1d^2;
    psn_max = power_noise.calc_psn(psn,power,chip,tsv,rho_m,mu_m,chip.temperature);
    psn.noise = psn_max;
    mismatch_norm = psn.noise/psn.noise_target;
    dispstr = sprintf('\tpsn_runs: %d\tNpads: %d\tpsn_target: %d\tpsn_max: %d\tmismatch_norm: %.3g',psn_iterations,psn.Npads, psn.noise_target, psn_max,mismatch_norm);
	disp(dispstr)
end

%% Final report

time_stop = cputime;
time_elapsed = time_stop - time_start;
disp(' ')
disp('Final block parameters:')
repstr = sprintf('\tNg_nom %d \t Ng_act: %d \t Atsv_nom: %.3g \t Atsv_act: %.3g \n\tN_tsvs: %d \t Npads_pow %d \t psn_nom %.4g \t psn_act %.4g', ...
                  chip.num_gates, chip.Ng_actual, tsv.max_area_fraction, tsv.actual_area_fraction, tsv.num, psn.Npads, psn.noise_target, psn_max);
disp(repstr)
repstr = sprintf('\th_tsv_um: %.4g \t w_tsv_um: %.4g',h_tsv_m/1e-6,w_tsv_m/1e-6);
disp(repstr);

disp(sprintf('Total time elapsed for block design: %d',time_elapsed))

