%the testbench for PDN; set the systimatical parameters
%this file lists all the parameters needed
%spparms('spumoni', 0);
chip.die = 100e-6; %die thickness
chip.bump = 25e-6; %bump thickness between two dice
%%%%%%%%%%%%%%%%%%for TSV domain%%%%%%%%%%%%%%%%%%%%%
TSV.d = 10e-6; %TSV diameter
TSV.px = 400e-6; %TSV X pitch
TSV.py = 400e-6; %TSV Y pitch

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

alpha = 0.00393;
TSV.rou = 17.1e-9*(1+alpha*(100-20)); % 17.1e-9 is the resistivity of copper @20 centigrade
TSV.R = TSV.rou*(chip.die+chip.bump)/(0.25*pi*TSV.d^2); %TSV resistance
TSV.L = 0.25*1.256e-6/pi*(chip.die+chip.bump)*log(TSV.px/TSV.d);
%%%%%%%%%%%%%%%%%for Chip parameters%%%%%%%%%%%%%%%%%
chip.Xsize = 1e-2; %chip x dimension; y dimension below
chip.Ysize = 1e-2;

chip.wire_thick = 5e-6; %wire thickness
chip.wire_p = 17.1e-9*(1+alpha*(100-20)); % wire resistivity
chip.wire_ar = 1.5; %aspect ratio = thickness/width (wire)
chip.wire_width = chip.wire_thick/chip.wire_ar;
chip.Rs = chip.wire_p/chip.wire_thick; %sheet resistance

chip.N = 2; %die number
chip.cap_per = [0.1 0.1];% 0.1 0.1]; %decap number
chip.cap_th = 0.9e-9; %capacitance effective thicknee (used for capacitance value calculation)
chip.c_gate = 3.9*8.85e-12/chip.cap_th;

chip.power = [2.82 74.49];% 74.49]; %74.49 74.49];  %total power dissipation of each die
%chip.power = [74.49 74.49]; %74.49 74.49];  %total power dissipation of each die
chip.map = [
            %%{
            0.2e-3 0.4e-3 4.1e-3 2.0e-3 0.4100 0.1; 
            0.2e-3 2.8e-3 4.1e-3 2.0e-3 0.3280 0.1;
            0.2e-3 5.2e-3 4.1e-3 2.0e-3 0.2460 0.1;
            0.2e-3 7.6e-3 4.1e-3 2.0e-3 0.1640 0.1;
            4.5e-3 0.2e-3 1.1e-3 9.7e-3 0.3735 0.1;
            5.8e-3 0.4e-3 4.1e-3 2.0e-3 0.2460 0.1;
            5.8e-3 2.8e-3 4.1e-3 2.0e-3 0.3690 0.1;
            5.8e-3 5.2e-3 4.1e-3 2.0e-3 0.1886 0.1;
            5.8e-3 7.6e-3 4.1e-3 2.0e-3 0.2542 0.1;
            %}
            0.2e-3 0.2e-3 1.0e-3 4.6e-3 2.3000 0.1;
            0.2e-3 5.0e-3 1.0e-3 3.6e-3 2.2690 0.1;
            1.4e-3 0.2e-3 7.3e-3 2.6e-3 14.8044 0.1;
            1.4e-3 3.0e-3 1.4e-3 5.6e-3 7.5264 0.1;
            3.0e-3 3.0e-3 1.4e-3 5.6e-3 8.3720 0.1;
            4.5e-3 3.0e-3 1.0e-3 3.5e-3 1.9600 0.1;
            4.5e-3 6.7e-3 1.0e-3 1.9e-3 1.3680 0.1;
            5.7e-3 3.0e-3 1.4e-3 5.6e-3 6.9776 0.1;
            7.3e-3 3.0e-3 1.4e-3 5.6e-3 10.5840 0.1;
            8.9e-3 0.2e-3 1.0e-3 4.2e-3 2.4360 0.1;
            8.9e-3 4.6e-3 1.0e-3 4.0e-3 3.0800 0.1;
            0.2e-3 8.8e-3 9.7e-3 1.0e-3 6.4020 0.1];
            %{
            0.2e-3 0.2e-3 1.0e-3 4.6e-3 2.3000 0.1;
            0.2e-3 5.0e-3 1.0e-3 3.6e-3 2.2690 0.1;
            1.4e-3 0.2e-3 7.3e-3 2.6e-3 14.8044 0.1;
            1.4e-3 3.0e-3 1.4e-3 5.6e-3 7.5264 0.1;
            3.0e-3 3.0e-3 1.4e-3 5.6e-3 8.3720 0.1;
            4.5e-3 3.0e-3 1.0e-3 3.5e-3 1.9600 0.1;
            4.5e-3 6.7e-3 1.0e-3 1.9e-3 1.3680 0.1;
            5.7e-3 3.0e-3 1.4e-3 5.6e-3 6.9776 0.1;
            7.3e-3 3.0e-3 1.4e-3 5.6e-3 10.5840 0.1;
            8.9e-3 0.2e-3 1.0e-3 4.2e-3 2.4360 0.1;
            8.9e-3 4.6e-3 1.0e-3 4.0e-3 3.0800 0.1;
            0.2e-3 8.8e-3 9.7e-3 1.0e-3 6.4020 0.1];
            %}
            %{
            0.2e-3 0.2e-3 1.0e-3 4.6e-3 2.3000 0.1;
            0.2e-3 5.0e-3 1.0e-3 3.6e-3 2.2690 0.1;
            1.4e-3 0.2e-3 7.3e-3 2.6e-3 14.8044 0.1;
            1.4e-3 3.0e-3 1.4e-3 5.6e-3 7.5264 0.2;
            3.0e-3 3.0e-3 1.4e-3 5.6e-3 8.3720 0.2;
            4.5e-3 3.0e-3 1.0e-3 3.5e-3 1.9600 0.1;
            4.5e-3 6.7e-3 1.0e-3 1.9e-3 1.3680 0.1;
            5.7e-3 3.0e-3 1.4e-3 5.6e-3 6.9776 0.2;
            7.3e-3 3.0e-3 1.4e-3 5.6e-3 10.5840 0.4;
            8.9e-3 0.2e-3 1.0e-3 4.2e-3 2.4360 0.1;
            8.9e-3 4.6e-3 1.0e-3 4.0e-3 3.0800 0.1;
            0.2e-3 8.8e-3 9.7e-3 1.0e-3 6.4020 0.1;          
            
            0.2e-3 0.2e-3 1.0e-3 4.6e-3 2.3000 0.1;
            0.2e-3 5.0e-3 1.0e-3 3.6e-3 2.2690 0.1;
            1.4e-3 0.2e-3 7.3e-3 2.6e-3 14.8044 0.1;
            1.4e-3 3.0e-3 1.4e-3 5.6e-3 7.5264 0.2;
            3.0e-3 3.0e-3 1.4e-3 5.6e-3 8.3720 0.2;
            4.5e-3 3.0e-3 1.0e-3 3.5e-3 1.9600 0.1;
            4.5e-3 6.7e-3 1.0e-3 1.9e-3 1.3680 0.1;
            5.7e-3 3.0e-3 1.4e-3 5.6e-3 6.9776 0.2;
            7.3e-3 3.0e-3 1.4e-3 5.6e-3 10.5840 0.4;
            8.9e-3 0.2e-3 1.0e-3 4.2e-3 2.4360 0.1;
            8.9e-3 4.6e-3 1.0e-3 4.0e-3 3.0800 0.1;
            0.2e-3 8.8e-3 9.7e-3 1.0e-3 6.4020 0.1];
            %}
%blk_num is for splitting the power maps of each die
% format: xl corner, yb corner, x size, y size, power, decap %
chip.map(:,6) = chip.map(:,6)*chip.c_gate;
%chip.map(13:24,5) = chip.map(13:24,5);

chip.blk_num = [9 12];% 12 12];

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
system.Vdd = 0.9; %VDD value for each pad

system.T = 50e-9; %transient simulation time 
system.Tr = 1e-9; %rise time
system.dt = 0.1e-9; %simulation step, TR scheme is stable for 0.1e-9, which is pretty optimized
system.ratio = 5; %the ratio between the maximum time step and minimum one
system.start_ratio = 1; %the start power of processor die; can optimize to generate more complicated power excitation
system.step = (1 - system.start_ratio)/(system.Tr/system.dt); %for each step, the increasement of the power excitation for processor die

power_noise_sim(chip, system, TSV);
%main simulation function
