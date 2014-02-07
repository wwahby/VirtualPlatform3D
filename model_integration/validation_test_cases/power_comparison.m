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
