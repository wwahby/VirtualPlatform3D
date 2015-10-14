%function plot_sweep_data( sweep, sweep_data, simulation) % just use as a
%script for now
fignum = 1;

%% EPC vs tiers

colors = [0 0 0; 0 0 1; 0 1 0; 1 0 0];
freq_mat = zeros(num_cooling_configs, num_stacks);
power_mat = zeros(num_cooling_configs, num_stacks);
for cind = 1:num_cooling_configs
    for nind = 1:num_stacks
        freq_mat(cind, nind) = freq(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
        power_mat(cind, nind) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
    end
end

Tclk_mat = 1./freq_mat;
epc_mat = power_mat .* Tclk_mat;
epc_mat = epc_mat'*1e9;

figure(1)
clf
hold on
b = bar(1:num_stacks, epc_mat, 1, 'grouped');
set(gca,'xtick',1:num_stacks)
set(gca,'xticklabel',tiers)
b(1).FaceColor = 'b';
b(2).FaceColor = 'y';
xlabel('Tiers')
ylabel('Energy Per Cycle (nJ)')
fixfigs(1,3,14,12)

figure(2)
clf
hold on
plot(1:num_stacks, epc_mat(:,1));
set(gca,'xtick',1:num_stacks)
set(gca,'xticklabel',tiers)
b(1).FaceColor = 'b';
xlabel('Tiers')
ylabel('Energy Per Cycle (nJ)')
fixfigs(2,3,14,12)