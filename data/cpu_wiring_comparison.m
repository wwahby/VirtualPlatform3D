%% Merom
intel_merom_pitch = [ 210 210 220 280 330 480 720 1080 ];
intel_merom_thickness = [ 170 190 200 250 300 430 650 975];
intel_merom_ar = [ 1.6 1.8 1.8 1.8 1.8 1.8 1.8 1.8];

deepak_merom_pitch = [ 220 220 283 283 283 283 880 880 ];
my_merom_n2_pitch = [ 220 220 410.5 410.5 740.9 740.9 ];
my_merom_n1_pitch = [ 220 220 267.4 541.6 740.9 ];
my_merom_n1_e03 = [ 220.0000	220.0000	220.0000	273.9328	483.7781	740.9495 ];
my_merom_n1_dyneff = [ 220.0000	220.0000	220.0000	267.3597	524.3668	746.8653 ];

my_merom_tdwithrep = [ 220.0000  220.0000  220.0000  230.1247  377.1137 ];
my_merom_tdwlarep2 = [ 220.0000  220.0000  220.0000  220.0000  255.7531  471.6866];
my_merom_sbbest = [0.2200    0.2200    0.2200    0.2850    0.4774    0.6771    0.8903    1.1446    1.4916] * 1e3;

% figure(1)
% clf
% plot(intel_merom_pitch,'k')
% hold on
% plot(deepak_merom_pitch,'r')
% plot(my_merom_n2_pitch,'b')
% plot(my_merom_n1_pitch,'g')
% plot(my_merom_n1_e03,'color',[0.7 0.3 0.7])
% plot(my_merom_n1_dyneff,'color',[0.3 0.7 1])
% xlabel('metal layer')
% ylabel('Pitch (nm)')
% fixfigs(1,3,14,12)

figure(7)
clf
plot(intel_merom_pitch,'k')
hold on
plot(deepak_merom_pitch,'g')
plot(my_merom_sbbest,'r')
ylim([0 1500])
xlabel('metal layer')
ylabel('Pitch (nm)')
title('Merom (65nm)')
fixfigs(7,3,14,12)

%% Penryn (45nm)
intel_penryn_pitch = [160 160 160 250 280 360 560 810 30500];
my_penryn_pitch = [160.0000  160.0000  160.0000  160.0000  160.0000  160.0000  172.1319  326.4828];
my_penryn_sbbest = [160.0000  160.0000  160.0000  160.0000  215.4635  328.8653  446.2047  572.0492  724.5117  933.4803];

figure(8)
clf
plot(intel_penryn_pitch,'k')
hold on
plot(my_penryn_sbbest,'r')
ylim([0 1000])
xlabel('metal layer')
ylabel('Pitch (nm)')
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

my_sb_tdwithrep = [ 112.5000  112.5000  112.5000  119.2457  232.3208  314.6503  397.2216];
my_sb_tdwlarep2 = [ 112.5000  112.5000  112.5000  112.5000  121.0211  232.5093 ];
my_sb_tdwlarep2best = [112.5000  112.5000  112.5000  112.5000  208.8135  311.7498  421.1808  550.4639  735.2589];

% figure(2)
% clf
% plot(intel_sb_pitch,'k')
% hold on
% plot(my_sb_n2_pitch,'b')
% plot(my_sb_n1_pitch,'g')
% plot(my_sb_n1_e03,'color',[0.7 0.3 0.7])
% plot(my_sb_n1_dyneff,'color',[0.3 0.7 1])
% xlabel('metal layer')
% ylabel('Pitch (nm)')
% ylim([0 1.1*max(my_sb_n2_pitch)])
% fixfigs(2,3,14,12)

figure(9)
clf
plot(intel_sb_pitch,'k')
hold on
plot(my_sb_tdwithrep,'b')
plot(my_sb_tdwlarep2best,'r')
xlabel('metal layer')
ylabel('Pitch (nm)')
ylim([0 1.1*max(my_sb_n2_pitch)])
title('Sandy Bridge (32nm)')
fixfigs(9,3,14,12)

%% Intel Ivy bridge

intel_ib_pitch =[	90	80	80	112	160	240	320	360	14000 ];
my_ib_n1_dyneff_pitch = [	90	90	90	115.836	210.5411	322.5414	493.471	];
my_ib_tdwithrep_pitch1 = [ 90.0000   90.0000   90.0000   90.0000   93.0672  164.8177  218.1233  249.1480 ];
my_ib_tdwithrep_pitch2 = [ 90.0000   90.0000   90.0000  112.9949  183.4550  249.1480];
my_ib_tdwlari3 = [90.0000   90.0000   90.0000   90.0000  123.0972  188.4797  257.9700  339.7838  457.6400];
my_ib_sbbest = [0.0900 0.0900 0.0900 0.0900 0.1610 0.2402 0.3196 0.3991 0.4794 0.5611 0.6463 0.7386 0.8447 0.9918 1.2331 ] * 1e3;

figure(10)
clf
plot(intel_ib_pitch,'k')
hold on
plot(my_ib_tdwithrep_pitch2,'b')
%plot(my_ib_tdwlari3,'r')
plot(my_ib_sbbest,'r')
xlabel('metal layer')
ylabel('Pitch (nm)')
ylim([0 1.1*max(my_ib_n1_dyneff_pitch)])
title('Ivy Bridge (22nm)')
fixfigs(10,3,14,12)

%% Num metal layers

intel_nodes = [65 45 32 22];
intel_nodes = [1 2 3 4];
intel_num_metal_layers = [8 9 8 9];
my_num_metal_layers = [6 8 6 9];
my_num_metal_layers_sbbest = [9 10 9 15];

figure(11)
clf
plot(intel_nodes,intel_num_metal_layers,'k')
hold on
plot(intel_nodes,my_num_metal_layers,'r')
plot(intel_nodes,my_num_metal_layers_sbbest,'m')
xlabel('Process node (nm)')
ylabel('Number of metal layers')
set(gca,'xtick',[1 2 3 4])
set(gca,'xticklabel',['65'; '45'; '32'; '22'])
ylim([5 15])
fixfigs(11,3,14,12)


%% Power
intel_nodes = [1 2 3 4];
intel_tdp = [44 55 95 77];
my_core_power = [21.8 24.7 21.74 69.86];
intel_num_cores = [2 2 4 4];
my_power_nomem = my_core_power*intel_num_cores;

figure(12)
clf
plot(intel_tdp,'k')
hold on
plot(intel_nodes,my_power_nomem,'r')
xlabel('Process node (nm)')
ylabel('Power (W)')
set(gca,'xtick',[1 2 3 4])
set(gca,'xticklabel',[65 45 32 22])


