function [max_temp temp_vec] = get_stack_temperature(num_layers,die_thickness,wire,tsv,chip_width,chip_height,package_width,package_height,heat,simulation,map,blk_num,power_therm_vec)

%% Thermal module -- Find actual system temperature

%%%%%%%%%%%%%%%%%%%%%%%geometry information of the chip%%%%%%%%%%%%%%%%%%%
    die.N = num_layers;
    die.model = 3; % for each die, how many layers we model
    if (die_thickness < heat.monolithic_max_thickness)
        die.material_IDs = heat.monolithic_material_IDs;
        die.layer_thicknesses = [die_thickness sum(wire.pn) heat.monolithic_intertier_bond_thickness];
        under_thickness = heat.monolithic_intertier_bond_thickness; % thickness of the bonding material between dice (underfill in this case)
    else
        die.material_IDs = heat.material_IDs; % Chip, ILD, underfill
        die.layer_thicknesses = [die_thickness sum(wire.pn) heat.underfill_thickness ];
        under_thickness =  heat.underfill_thickness; % thickness of the bonding material between dice (oxide in this case)
    end
    

    % flip chip package; order:
    %heatsink->TIM->CHIP_BULK1->METAL->BONDING->CHIP_BULK2-> ...
    %               CHIP_BULKN->METAL->MICRO-BUMPS->INTERPOSER    
    thick.bump = heat.bump_thickness;  %micro-bump thickness; between second die and interposer
    thick.tim = heat.tim_thickness; %tim thickness; between the chip and heatsink
    thick.under = under_thickness; %underfill bonding thickness; between two dies
    thick.inter = heat.interposer_thickness; %interposer thickness
    
    thick.die = die_thickness; %die thickness
    thick.ild = sum(wire.pn); %metal layer thickness

    grid_factor = 50;
    chip_therm.Xsize = chip_width; %x dimension of chip
    chip_therm.Ysize = chip_height; %y dimension of chip
    chip_therm.Xgrid = chip_therm.Xsize/grid_factor; %x grid size of chip
    chip_therm.Ygrid = chip_therm.Ysize/grid_factor; %y grid size of chip
    
    %for the interposer dimension
    pack.Xsize = package_width;
    pack.Ysize = package_height;
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
    
    %power_therm_vec = ones(1,chip.num_layers)*chip_power/chip.num_layers;  %power dissipation of each die
    % from top to bottom; unit: watt
    
    granularity = 20;
    % thermal map, number of color used
    
    draw = simulation.draw_thermal_map;
    drawP = simulation.draw_thermal_map;
    % whether to draw the thermal map; 1 yes; 0 no
    
    displayT = simulation.print_thermal_data;
    % print the temperature information


    temp_vec = thermal.ThermSim( die, thick, chip_therm, pack, ...
              tsv_therm, bump, portion, power_therm_vec, ...
              map, blk_num, granularity, draw, drawP, heat, displayT);
          
    max_temp = max(temp_vec);