%function plot_sweep_data( sweep, sweep_data, simulation) % just use as a
%script for now
fignum = 1;

%% EPC, Freq, Temp, Total Power, Total Power Density

colors = [0 0 0; 0 0 1; 0 1 0; 1 0 0];
freq_mat = zeros(num_stacks, num_cooling_configs);
power_mat = zeros(num_stacks, num_cooling_configs);
power_density_mat = zeros(num_stacks, num_cooling_configs);
wire_power_mat = zeros(num_stacks, num_cooling_configs);
rep_power_mat = zeros(num_stacks, num_cooling_configs);
dynamic_power_mat = zeros(num_stacks, num_cooling_configs);
leakage_power_mat = zeros(num_stacks, num_cooling_configs);
area_mat = zeros(num_stacks, num_cooling_configs);

for nind = 1:num_stacks
    for cind = 1:num_cooling_configs
        freq_mat(nind, cind) = freq(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
        power_mat(nind, cind) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
        power_density_mat(nind, cind) = power_density(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
        wire_power_mat(nind, cind) = wire_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
        rep_power_mat(nind, cind) = rep_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
        dynamic_power_mat(nind, cind) = dynamic_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
        leakage_power_mat(nind, cind) = leakage_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
        area_mat(nind, cind) = chip_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind}.area_per_layer_m2;                               
    end
end
Tclk_mat = 1./freq_mat;
epc_mat = power_mat .* Tclk_mat;
epc_mat = epc_mat'*1e9;



%%

figure(1)
clf
b = bar(freq_mat/1e9, 1, 'grouped');
b(1).FaceColor = 'b';
b(2).FaceColor = 'y';
xlabel('Tiers')
ylabel('Maximum Frequency (GHz)')
fixfigs(1,3,14,12)