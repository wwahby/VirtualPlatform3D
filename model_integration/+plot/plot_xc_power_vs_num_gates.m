
power_tot_vec = zeros(1, num_gate_sweeps);
power_logic_vec = zeros(1, num_gate_sweeps);
power_xc_vec = zeros(1, num_gate_sweeps);
area_tot_vec = zeros(1,num_gate_sweeps);
area_per_layer_vec = zeros(1,num_gate_sweeps);

power_logic_mat = logic_dynamic_power + logic_leakage_power;         

for num_gates_ind = 1:num_gate_sweeps
    power_tot_vec(num_gates_ind) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind,num_gates_ind);
    power_logic_vec(num_gates_ind) = power_logic_mat(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind,num_gates_ind);
    area_per_layer_vec(num_gates_ind) = area_per_layer_mat(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind,num_gates_ind);
    area_tot_vec(num_gates_ind) = area_tot_mat(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind,num_gates_ind);
    power_xc_vec(num_gates_ind) = interconnect_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind,num_gates_ind);
    
end
power_xc_vec = power_tot_vec - power_logic_vec;
power_density_xc_vec = power_xc_vec ./ area_per_layer_vec;
power_density_tot_vec = power_tot_vec ./ area_per_layer_vec;
power_density_logic_vec = power_logic_vec ./ area_per_layer_vec;

colors = {'k', 'b', 'g', 'r'};

figure(1)
clf
hold on
plot(num_gates_vec, power_tot_vec,'k')
plot(num_gates_vec, power_xc_vec,'b')
%plot(num_gates_vec, power_logic_vec, 'r')
xlabel('Number of gates')
ylabel('Power (W)')
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'xtick', 10.^(1:2:9))
set(gca,'ytick', 10.^(-6:2:4))
xlim([1e1 1e9])
grid on
fixfigs(1,3,14,12)

figure(2)
clf
hold on
plot(num_gates_vec, power_density_tot_vec, 'k')
plot(num_gates_vec, power_density_xc_vec, 'b')
%plot(num_gates_vec, power_density_logic_vec, 'r')
xlabel('Number of gates')
ylabel('Power density (W/m^2)')
set(gca,'xscale','log')
set(gca,'yscale','log')
set(gca,'xtick', 10.^(1:1:9))
xlim([1e2 1e9])
grid on
fixfigs(2,3,14,12)

