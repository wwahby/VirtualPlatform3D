%name = {	backprop	fastwalsh	gaussian	heartwall	matrixmul	particlefilter	similarityscore	s.reduce	s.srad2	stringmatch	geomean };
b0p5 = [	0.78	0.61	0.92	0.98	0.98	0.88	0.99	0.99	0.76	0.98	0.89];
b2p0 = [	1.04	1.5	1.08	1.14	1.44	1.24	1.48	1.02	1.78	1.2	1.25 ];
b4p0 = [	1.05	1.74	1.1	1.19	1.82	1.38	1.79	1.04	2	1.28	1.38 ];
b8p0 = [	1.05	1.76	1.12	1.21	1.92	2	1.88	1.04	2.02	1.32	1.4 ];
binf = [	1.15	1.96	1.9	1.24	13.7	2.5	13.5	1.39	2.6	4	2.6 ];



figure(1)
clf
hold on
grid on
% cdfplot(b0p5)
% cdfplot(b2p0)
% cdfplot(b4p0)
% cdfplot(b8p0)
cdfplot(binf)
fixfigs(1,2,14,12)
% set(gca,'xscale','log')
% xlim([5e-1 2e1])
xlim([1e0, 5])
ylim([0, 1])
xlabel('Speedup over baseline (X)')
ylabel('Fraction of benchmarks')
title('')
set(gca,'ytick', 0:0.1:1);


%%

tiers_act = [ 8 8 8 8 8 ];
tiers_pred = [8 8 7 6 8];

power_act = [ 65 65 95 73 95];
power_pred = [60.27 63.62 105.56 52.74 91.80 ];

names = {'65nm B', '45nm A', '45nm B', '32nm A', '32nm B'};

figure(2)
clf
hold on
grid on
plot(tiers_act, 'b')
plot(tiers_pred, 'r')
ylabel('Signal wiring tiers')
set(gca, 'xtick', 1:length(names));
set(gca, 'xticklabel', names);
fixfigs(2,2,14,12)

figure(3)
clf
hold on
grid on

b = bar( [ power_act ; power_pred ]', 1, 'grouped');
ylabel('Power (W)')
set(gca, 'xtick', 1:length(names));
set(gca, 'xticklabel', names)
xlim([0.5, length(names)+0.5])
b(1).FaceColor = 'b';
b(2).FaceColor = 'y';
fixfigs(3,2,14,12)