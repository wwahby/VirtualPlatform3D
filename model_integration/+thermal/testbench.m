%%%%%%%%%%%%%%%%%%%%%%%%geometry information of the chip%%%%%%%%%%%%%%%%%%%
    die.N = 4;
    die.model = 3; % for each die, how many layers we model
    
    % flip chip package; order:
    %heatsink->TIM->CHIP_BULK1->METAL->BONDING->CHIP_BULK2-> ...
    %               CHIP_BULKN->METAL->MICRO-BUMPS->INTERPOSER    
    thick.bump = 40e-6;  %micro-bump thickness; between second die and interposer
    thick.tim = 5e-6; %tim thickness; between the chip and heatsink
    thick.under = 5e-6; %underfill bonding thickness; between two dies
    thick.inter = 200e-6; %interposer thickness
    thick.die = 50e-6; %die thickness
    thick.ild = 5e-6; %metal layer thickness
    
    
    chip.Xsize = 10e-3; %x dimension of chip
    chip.Ysize = 10e-3; %y dimension of chip
    chip.Xgrid = 200e-6; %x grid size of chip
    chip.Ygrid = 200e-6; %y grid size of chip
    
    %for the interposer dimension
    pack.Xsize = 35e-3;
    pack.Ysize = 35e-3;
    pack.Xgrid = 200e-6;
    pack.Ygrid = 200e-6;
    
    %assumed TSV starting from top metal layer of a top die
    %            to the first metal layer of a bottom die
    %Thus the TSV passes through bonding layer & bulk of a botom die
    %TSV geometry 
    tsv.d = 5e-6; % tsv diameter including the liner thickness
    tsv.liner = 0.5e-6; %liner thickness
    tsv.px = 50e-6; %x direction pitch
    tsv.py = 50e-6; %y direction pitch
    tsv.Nx = 200; %x direction number
    tsv.Ny = 200; %y direction number

    %Bump geometry
    bump.d = 20e-6;
    bump.px = 100e-6;
    bump.py = 100e-6;
    bump.Nx = 40;
    bump.Ny = 40;
    
    portion = 0.5; %the metal portion in the ILD layers
    %This is used for equivalent thermal resistance calculation of metal
    %layer
    power = [40, 30, 20, 10];  %power dissipation of each die
    % from top to bottom; unit: watt
    
    %power blocks by each die
    %format: bottom left point bl_x, bl_y, width, height, power
    %list blocks in die1 and then die2, die3 ....
    map = [0     0     5e-3    5e-3     8;
           0     5e-3  5e-3    5e-3     2;
           5e-3  0     5e-3    10e-3    30;
           0     0     8e-3    6e-3     20;
           0     6e-3  8e-3    4e-3     4;
           8e-3  0     2e-3    10e-3    6;
           1e-3  1e-3  8e-3    8e-3     13];
    %blk_num is for splitting the power maps of each die
    blk_num = [3 3 1 0];
    
    
    granularity = 20;
    % thermal map, number of color used
    draw = 1;
    % whether to draw the thermal map; 1 yes; 0 no
    drawP = 0;
    % whether to draw the power maps; 1 yes, 0 no
    % only draw the die with blk_num > 1 (non-uniform cases)
    displayT = 1;
    % print the temperature information
%%%%%%%%%%%%%%%%%%%%%%%%finish geometry information%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%boundary condition%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %the thermal resistivity of each boundary
    % r = 1/(hA); A is the size of top surface area
    % the cooling capability of the top heatsink; 20000, 1cm*1cm, means:
    % 0.5 W*cm^2/K
    h.up = 20000;
    
    h.down = 0; % the cooling of bottom surface; 
    %(only the area with the same size of chip;)
    %microfluidic is assumed to be as large as chip in the interposer
    
    h.side = 0;
    % side surface cooling, usually near adiabatic
    
    h.d = 0;
    %the cooling of the bottom surface except for the MFHS area
    
    h.Ta = 298;
    %the ambient temperature
%%%%%%%%%%%%%%%%%%%%%%%%%finish boundary condition%%%%%%%%%%%%%%%%%%%%%%%%%
    
    T = thermal.ThermSim( die, thick, chip, pack, ...
              tsv, bump, portion, power,  ...
              map, blk_num, granularity, ...
              draw, drawP, h, displayT);
               
               