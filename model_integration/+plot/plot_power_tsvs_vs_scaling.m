%function plot_sweep_data( sweep, sweep_data, simulation) % just use as a
%script for now
fignum = 1;

%% EPC, Freq, Temp, Total Power, Total Power Density

colors = [0 0 0; 0 0 1; 0 1 0; 1 0 0];
num_metal_levels_mat = zeros(num_wire_resistivities, num_scaling_factors);
npads_mat = zeros(num_wire_resistivities, num_scaling_factors);

for nind = 1
    for wire_res_ind = 1:num_wire_resistivities
        for scaling_ind = 1:num_scaling_factors
            if (routable_design(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind) == 1)
                
                npads_mat(wire_res_ind, scaling_ind) = npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
                num_metal_levels_mat(wire_res_ind, scaling_ind) = num_metal_levels(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
            else
                npads_mat(wire_res_ind, scaling_ind) = 0;
                num_metal_levels_mat(wire_res_ind, scaling_ind) = 0;
            end
                
        end
    end
end

figure(1)
clf
hold on
b = bar(npads_mat', 1, 'grouped');
b(1).FaceColor = 'b';
b(2).FaceColor = 'y';
% b(3).FaceColor = 'y';
% b(4).FaceColor = 'r';
xlim([0.5, num_scaling_factors + .5])
set(gca, 'xtick', 1:num_scaling_factors)
set(gca, 'xticklabel', node_labels)
xlabel('Process Node')
ylabel('Number of Power TSVs/Pads')
fixfigs(1,3,14,12)

%set(gca,'yscale','log')



