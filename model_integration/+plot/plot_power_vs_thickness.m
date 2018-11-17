
colors = {'k', 'b', 'r', 'g'};

figure(1)
clf
hold on
grid on
for nind = 1:num_stacks
    power_vec = zeros(1,num_thicks);
    power_density_vec = zeros(1,num_thicks);
    xc_power_vec = zeros(1,num_thicks);
    for thind=1:num_thicks
        power_vec(thind) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
        power_density_vec(thind) = power_density(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);   
        xc_power_vec(thind) = interconnect_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind,num_gates_ind);
        
    end
    plvec = power_vec;
        
    plot(thicknesses*1e6, plvec, colors{nind}, 'linewidth', 1.5)
end

set(gca,'xscale','log')
xlabel('Tier Thickness (microns)')
ylabel('Power (W)')