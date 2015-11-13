

num_metal_levels_vec = zeros(1,num_perms);
power_vec = zeros(1,num_perms);
for pind = 1:num_perms
    num_metal_levels_vec(pind) = num_metal_levels(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
    power_vec(pind) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);                                                 
end

figure(1)
clf
hold on
plot(rel_permittivities, num_metal_levels_vec,'k');
plot(rel_permittivities, power_vec, 'r');
