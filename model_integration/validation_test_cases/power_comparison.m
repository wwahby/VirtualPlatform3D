% simulated vs. actual full-die power consumption

intel_node = [ 65 45 45 32 ];
intel_node_plot = [1 2 3 4];
intel_tdp = [65 65 95 73];
actual_power = [68 57.5 94.6 67.8];

figure(1)
clf
plot(intel_node_plot,intel_tdp,'k')
hold on
plot(intel_node_plot,actual_power,'r')
set(gca,'xtick',[1 2 3 4])
set(gca,'xticklabel',[65 45 45 32])
set(gca,'xticklabel',{ '65nm Core' '45nm Core' '45nm Nehalem' '32nm Nehalem' })
fixfigs(1,3,14,12)



%% Updated comparison
% simulated vs. actual full-die power consumption

intel_node = [ 65 45 45 32 32 ];
intel_node_plot = [1 2 3 4 5];
intel_tdp = [65 65 95 73 95 ];
sim_power = [58.88 59.49 103.11 59.18 116];

intel_tmax_core = [100 105 100 105 105];
sim_temp = [74.31 81.94 113.5 68.85 90];
%sim_temp = [74.31 81.94 60.57 68.85 90];

figure(2)
clf
plot(intel_node_plot,intel_tdp,'k')
hold on
plot(intel_node_plot,sim_power,'r')
set(gca,'xtick',[1 2 3 4 5])
set(gca,'xticklabel',[65 45 45 32 32])
set(gca,'xticklabel',{ '65nm Core' '45nm Core' '45nm Nehalem' '32nm Nehalem' '32nm SB'})
ylabel('Power (W)')
fixfigs(2,3,14,12)

figure(3)
clf
plot(intel_node_plot,intel_tmax_core,'k')
hold on
plot(intel_node_plot,sim_temp,'r')
set(gca,'xtick',[1 2 3 4 5])
set(gca,'xticklabel',[65 45 45 32 32])
set(gca,'xticklabel',{ '65nm Core' '45nm Core' '45nm Nehalem' '32nm Nehalem' '32nm SB'})
fixfigs(3,3,14,12)