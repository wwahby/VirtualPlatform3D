

num_metal_levels_vec = zeros(num_stacks,num_perms);
power_vec = zeros(num_stacks,num_perms);
for nind = 1:num_stacks
    for pind = 1:num_perms
        num_metal_levels_vec(nind,pind) = num_metal_levels(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
        power_vec(nind, pind) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);                                                 
    end
end

colors = {'k', 'b', 'g', 'r'};

figure(1)
clf
hold on
for nind = 1:num_stacks
    plot(rel_permittivities, num_metal_levels_vec(nind,:),'color', colors{nind}, 'linestyle', '-');
end
xlabel('Relative Permittivity')
ylabel('Number of Metal Levels')
fixfigs(1,3,14,12)

figure(2)
clf
hold on
for nind = 1:num_stacks
    plot(rel_permittivities, power_vec(nind,:), 'color', colors{nind}, 'linestyle', '-');
end
xlabel('Relative Permittivity')
ylabel('Power (W)')
fixfigs(2,3,14,12)