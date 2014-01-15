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
plot(my_merom_n1_dyneff,'b')
plot(my_merom_tdwithrep,'r')
xlabel('metal layer')
ylabel('Pitch (nm)')
fixfigs(7,3,14,12)

%% Sandy bridge
intel_sb_pitch = [ 112.5	112.5	112.5	168.8	225	337.6	450.1	566.5	19400 ];
intel_sb_thickness = [ 95 95 95 151 204 303 388 504 8000 2.55e4 ];
intel_sb_ar = [ 1.7 1.7 1.7 1.8 1.8 1.8 1.7 1.8 1.5 ];
my_sb_n2_pitch = [ 112.5	112.5	199.3	199.3	420	420	718.2	718.2 ];
my_sb_n1_pitch = [ 112.5	112.5	131.4947	247.3285	374.1441	540.3985	718.2076 ];
my_sb_n1_e03 = [ 112.5000	112.5000	112.5000	155.0817	252.8671	357.0507	480.4285	712.0960	712.7645 ];
my_sb_n1_dyneff = [ 112.5000	112.5000	112.5000	148.2061	265.4723	394.7707	577.4500	718.2076 ];

my_sb_tdwithrep = [112.5000  112.5000  112.5000  119.2457  232.3208  314.6503  397.2216];

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

figure(8)
clf
plot(intel_sb_pitch,'k')
hold on
plot(my_sb_n1_dyneff,'b')
plot(my_sb_tdwithrep,'r')
xlabel('metal layer')
ylabel('Pitch (nm)')
ylim([0 1.1*max(my_sb_n2_pitch)])
fixfigs(8,3,14,12)

%% Intel Ivy bridge

intel_ib_pitch =[	90	80	80	112	160	240	320	360	14000 ];
my_ib_n1_dyneff_pitch = [	90	90	90	115.836	210.5411	322.5414	493.471	];
my_ib_tdwithrep_pitch = [ 90.0000   90.0000   90.0000   90.0000   93.0672  164.8177  218.1233  249.1480 ];
my_ib_tdwithrep_pitch = [ 90.0000   90.0000   90.0000  112.9949  183.4550  249.1480];

figure(9)
clf
plot(intel_ib_pitch,'k')
hold on
plot(my_ib_n1_dyneff_pitch,'b')
plot(my_ib_tdwithrep_pitch,'r')
xlabel('metal layer')
ylabel('Pitch (nm)')
ylim([0 1.1*max(my_ib_n1_dyneff_pitch)])
fixfigs(9,3,14,12)

