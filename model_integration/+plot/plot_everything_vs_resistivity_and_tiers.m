

num_metal_levels_vec = zeros(num_stacks,num_wire_resistivities);
power_vec = zeros(num_stacks,num_wire_resistivities);
npads_vec = zeros(num_stacks,num_wire_resistivities);
for nind = 1:num_stacks
    for wire_res_ind = 1:num_wire_resistivities
        num_metal_levels_vec(nind,wire_res_ind) = num_metal_levels(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
        power_vec(nind, wire_res_ind) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
        npads_vec(nind, wire_res_ind) = npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
                                                    
    end
end

colors = {'k', 'b', 'g', 'r'};

figure(1)
clf
hold on
for nind = 1:num_stacks
    plot(wire_resistivities*1e9, num_metal_levels_vec(nind,:),'color', colors{nind}, 'linestyle', '-');
end
xlim([10 60])
xlabel('Bulk Resistivity (\Omeganm)')
ylabel('Number of Metal Levels')
xlim([10 60])
fixfigs(1,3,14,12)

figure(2)
clf
hold on
for nind = 1:num_stacks
    plot(wire_resistivities*1e9, power_vec(nind,:), 'color', colors{nind}, 'linestyle', '-');
end
xlabel('Bulk Resistivity (\Omeganm)')
ylabel('Power (W)')
xlim([10 60])
fixfigs(2,3,14,12)

figure(3)
clf
hold on
for nind = 1:num_stacks
    plot(wire_resistivities*1e9, npads_vec(nind,:), 'color', colors{nind}, 'linestyle', '-');
end
xlabel('Bulk Resistivity (\Omeganm)')
ylabel('Number of Power TSVs Per Die')
set(gca,'yscale','log')
xlim([10 60])
fixfigs(3,3,14,12)


total_die_area_m2 = core.chip.area_total;
total_die_area_mm2 = total_die_area_m2 * 1e6;
tsv_max_area_fraction = 0.01;

psn_tsv_area_vec_m2 = npads_vec * core.psn.power_tsv_width^2;
psn_tsv_area_vec_mm2 = psn_tsv_area_vec_m2*1e6;
figure(5)
clf
hold on
for nind = 1:num_stacks
    area_per_die_mm2 = total_die_area_mm2/tiers(nind);
    total_tsv_area_allocation_mm2 = area_per_die_mm2 * tsv_max_area_fraction;
    
    plot(wire_resistivities*1e9, psn_tsv_area_vec_mm2(nind,:)/total_tsv_area_allocation_mm2, 'color', colors{nind}, 'linestyle', '-');
end
xlabel('Bulk Resistivity (\Omeganm)')
ylabel('Power TSV Area / Max TSV Area')
xlim([10 60])
set(gca,'yscale','log')
fixfigs(5,3,14,12)