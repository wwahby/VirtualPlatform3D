%% Merom
intel_merom_pitch = [ 210 210 220 280 330 480 720 1080 ];
intel_merom_thickness = [ 170 190 200 250 300 430 650 975];
intel_merom_ar = [ 1.6 1.8 1.8 1.8 1.8 1.8 1.8 1.8];

deepak_merom_pitch = [ 220 220 283 283 283 283 880 880 ];
my_merom_pitch = [0.2200    0.2200    0.2200    0.2200    0.3753    0.5540    0.7602    1.0547]*1e3;

figure(7)
clf
plot(intel_merom_pitch,'k')
hold on
plot(deepak_merom_pitch,'g')
plot(my_merom_pitch,'r')
xlabel('metal layer')
ylabel('Pitch (nm)')
title('Merom (65nm)')
fixfigs(7,3,14,12)

%% Intel Penryn

intel_penryn_pitch = [160	160	160	250	280	360	560	810	30500];
my_penryn_pitch = [160.0000  160.0000  160.0000  160.0000  160.0000  228.1994  340.6842  470.0975  660.0702];

figure(8)
clf
plot(intel_penryn_pitch,'k')
hold on
plot(my_penryn_pitch,'r')
xlabel('metal layer')
ylabel('Pitch (nm)')
ylim([0 1.1*max(intel_penryn_pitch(1:end-1))])
title('Penryn (45nm)')
fixfigs(8,3,14,12)


%% Sandy bridge
intel_sb_pitch = [ 112.5	112.5	112.5	168.8	225	337.6	450.1	566.5	19400 ];
intel_sb_thickness = [ 95 95 95 151 204 303 388 504 8000 2.55e4 ];
intel_sb_ar = [ 1.7 1.7 1.7 1.8 1.8 1.8 1.7 1.8 1.5 ];
my_sb_n2_pitch = [ 112.5	112.5	199.3	199.3	420	420	718.2	718.2 ];
my_sb_n1_pitch = [ 112.5	112.5	131.4947	247.3285	374.1441	540.3985	718.2076 ];
my_sb_n1_e03 = [ 112.5000	112.5000	112.5000	155.0817	252.8671	357.0507	480.4285	712.0960	712.7645 ];
my_sb_n1_dyneff = [ 112.5000	112.5000	112.5000	148.2061	265.4723	394.7707	577.4500	718.2076 ];

my_sb_tdwithrep = [112.5000  112.5000  112.5000  119.2457  232.3208  314.6503  397.2216];
my_sb_pitch = [ 112.5000  112.5000  112.5000  145.7779  251.6180  368.1610  532.7681  657.6356];

figure(9)
clf
plot(intel_sb_pitch,'k')
hold on
plot(my_sb_pitch,'r')
xlabel('metal layer')
ylabel('Pitch (nm)')
ylim([0 1.1*max(intel_sb_pitch(1:end-1))])
title('Sandy Bridge (32nm)')
fixfigs(9,3,14,12)

%% Intel Ivy bridge

intel_ib_pitch =[	90	80	80	112	160	240	320	360	14000 ];
my_ib_pitch = [ 90.0000   90.0000   90.0000   98.4239  172.5822  254.3545  369.3061  457.6400];

figure(10)
clf
plot(intel_ib_pitch,'k')
hold on
plot(my_ib_pitch,'r')
xlabel('metal layer')
ylabel('Pitch (nm)')
ylim([0 1.1*max(intel_ib_pitch(1:end-1))])
title('Ivy Bridge (22nm)')
fixfigs(10,3,14,12)

%% Num metal layers

intel_nodes = [65 45 32 22];
intel_nodes = [1 2 3 4];
intel_num_metal_layers = [8 9 8 9];
my_num_metal_layers = [length(my_merom_pitch) length(my_penryn_pitch) length(my_sb_pitch) length(my_ib_pitch) ];

figure(11)
clf
plot(intel_nodes,intel_num_metal_layers,'k')
hold on
plot(intel_nodes,my_num_metal_layers,'r')
%plot(intel_nodes,my_num_metal_layers_sbbest,'m')
xlabel('Process node (nm)')
ylabel('Number of metal layers')
set(gca,'xtick',[1 2 3 4])
set(gca,'xticklabel',['65'; '45'; '32'; '22'])
ylim([5 10])
fixfigs(11,3,14,12)


%% Power
intel_nodes = [1 2 3 4];
intel_tdp = [44 55 95 77];
my_core_power = [21.8 24.7 21.74 69.86];
my_core_power = [31.5 46.98 45.11 71.7];
intel_num_cores = [2 2 4 4];
my_power_nomem = my_core_power.*intel_num_cores;

figure(12)
clf
plot(intel_tdp,'k')
hold on
plot(intel_nodes,my_power_nomem,'r')
xlabel('Process node (nm)')
ylabel('Power (W)')
set(gca,'xtick',[1 2 3 4])
set(gca,'xticklabel',[65 45 32 22])
fixfigs(12,3,14,12)


