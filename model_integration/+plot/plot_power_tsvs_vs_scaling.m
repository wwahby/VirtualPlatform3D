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

figure(2)
clf
hold on
b = bar(num_metal_levels_mat', 1, 'grouped');
b(1).FaceColor = 'b';
b(2).FaceColor = 'y';
% b(3).FaceColor = 'y';
% b(4).FaceColor = 'r';
xlim([0.5, num_scaling_factors + .5])
set(gca, 'xtick', 1:num_scaling_factors)
set(gca, 'xticklabel', node_labels)
xlabel('Process Node')
ylabel('Number of Metal Levels')
fixfigs(2,3,14,12)


%% ALL tiers

num_metal_levels_mat = zeros(num_wire_resistivities*num_stacks, num_scaling_factors);
npads_mat = zeros(num_wire_resistivities*num_stacks, num_scaling_factors);
power_mat = zeros(num_wire_resistivities*num_stacks, num_scaling_factors);
temp_mat = zeros(num_wire_resistivities*num_stacks, num_scaling_factors);
for nind = 1:num_stacks
    for wire_res_ind = 1:num_wire_resistivities
        for scaling_ind = 1:num_scaling_factors
            if (routable_design(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind) == 1)
                power_mat(wire_res_ind*num_stacks + (nind-1), scaling_ind) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
                temp_mat(wire_res_ind*num_stacks + (nind-1), scaling_ind) = temperature(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
                npads_mat(wire_res_ind*num_stacks + (nind-1), scaling_ind) = npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
                num_metal_levels_mat(wire_res_ind*num_stacks + (nind-1), scaling_ind) = num_metal_levels(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
            else
                power_mat(wire_res_ind*num_stacks + (nind-1), scaling_ind) = 0;
                temp_mat(wire_res_ind*num_stacks + (nind-1), scaling_ind) = 0;
                npads_mat(wire_res_ind*num_stacks + (nind-1), scaling_ind) = 0;
                num_metal_levels_mat(wire_res_ind*num_stacks + (nind-1), scaling_ind) = 0;
            end
                
        end
    end
end

figure(3)
clf
hold on
b = bar(num_metal_levels_mat', 1, 'grouped');
% b(3).FaceColor = 'y';
% b(4).FaceColor = 'r';
xlim([0.5, num_scaling_factors + .5])
set(gca, 'xtick', 1:num_scaling_factors)
set(gca, 'xticklabel', node_labels)
colormap jet
xlabel('Process Node')
ylabel('Number of Metal Levels')
fixfigs(3,3,14,12)

figure(4)
clf
hold on
b = bar(power_mat', 1, 'grouped');
% b(3).FaceColor = 'y';
% b(4).FaceColor = 'r';
xlim([0.5, num_scaling_factors + .5])
set(gca, 'xtick', 1:num_scaling_factors)
set(gca, 'xticklabel', node_labels)
colormap jet
xlabel('Process Node')
ylabel('Power (W)')
fixfigs(4,3,14,12)

figure(5)
clf
hold on
b = bar(temp_mat', 1, 'grouped');
% b(3).FaceColor = 'y';
% b(4).FaceColor = 'r';
xlim([0.5, num_scaling_factors + .5])
set(gca, 'xtick', 1:num_scaling_factors)
set(gca, 'xticklabel', node_labels)
colormap jet
xlabel('Process Node')
ylabel('Maximum Temperature (C)')
fixfigs(5,3,14,12)

figure(6)
clf
hold on
b = bar(npads_mat', 1, 'grouped');
xlim([0.5, num_scaling_factors + .5])
set(gca, 'xtick', 1:num_scaling_factors)
set(gca, 'xticklabel', node_labels)
colormap jet
xlabel('Process Node')
ylabel('Number of Power TSVs/Pads')
%set(gca,'yscale','log')
fixfigs(6,3,14,12)

%%

num_metal_levels_mat = zeros(num_wire_resistivities, num_stacks, num_scaling_factors);
npads_mat = zeros(num_wire_resistivities, num_stacks, num_scaling_factors);
power_mat = zeros(num_wire_resistivities, num_stacks, num_scaling_factors);
temp_mat = zeros(num_wire_resistivities, num_stacks, num_scaling_factors);
for nind = 1:num_stacks
    for wire_res_ind = 1:num_wire_resistivities
        for scaling_ind = 1:num_scaling_factors
            if (routable_design(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind) == 1)
                power_mat(wire_res_ind, nind, scaling_ind) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
                temp_mat(wire_res_ind, nind, scaling_ind) = temperature(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
                npads_mat(wire_res_ind, nind, scaling_ind) = npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
                num_metal_levels_mat(wire_res_ind, nind, scaling_ind) = num_metal_levels(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
            else
                power_mat(wire_res_ind, nind, scaling_ind) = 0;
                temp_mat(wire_res_ind, nind, scaling_ind) = 0;
                npads_mat(wire_res_ind, nind, scaling_ind) = 0;
                num_metal_levels_mat(wire_res_ind, nind, scaling_ind) = 0;
            end
                
        end
    end
end


colors = {'k', 'b', 'r', 'y'};
linestyles = {'-', '--', ':', '-.'};
figure(7)
clf
hold on
for wire_res_ind = 1:num_wire_resistivities
    for nind = 1:num_stacks
        plot_vec = zeros(1,num_scaling_factors);
        plot_vec(1,:) = num_metal_levels_mat(wire_res_ind,nind,:);
        plot_vec = plot_vec( plot_vec > 0);
        plot(1:length(plot_vec), plot_vec, 'color', colors{wire_res_ind}, 'linestyle', linestyles{nind})
    end
end
set(gca, 'xtick', 1:num_scaling_factors)
set(gca, 'xticklabel', node_labels)
xlabel('Process Node')
ylabel('Number of Metal Levels')
fixfigs(7,3,14,12)

figure(8)
clf
hold on
for wire_res_ind = 1:num_wire_resistivities
    for nind = 1:num_stacks
        plot_vec = zeros(1,num_scaling_factors);
        plot_vec(1,:) = power_mat(wire_res_ind,nind,:);
        plot_vec = plot_vec( plot_vec > 0);
        plot(1:length(plot_vec), plot_vec, 'color', colors{wire_res_ind}, 'linestyle', linestyles{nind})
    end
end
set(gca, 'xtick', 1:num_scaling_factors)
set(gca, 'xticklabel', node_labels)
colormap jet
xlabel('Process Node')
ylabel('Power (W)')
fixfigs(8,3,14,12)

figure(9)
clf
hold on
for wire_res_ind = 1:num_wire_resistivities
    for nind = 1:num_stacks
        plot_vec = zeros(1,num_scaling_factors);
        plot_vec(1,:) = temp_mat(wire_res_ind,nind,:);
        plot_vec = plot_vec( plot_vec > 0);
        plot(1:length(plot_vec), plot_vec, 'color', colors{wire_res_ind}, 'linestyle', linestyles{nind})
    end
end
set(gca, 'xtick', 1:num_scaling_factors)
set(gca, 'xticklabel', node_labels)
colormap jet
xlabel('Process Node')
ylabel('Maximum Temperature (C)')
fixfigs(9,3,14,12)

figure(10)
clf
hold on
for wire_res_ind = 1:num_wire_resistivities
    for nind = 1:num_stacks
        plot_vec = zeros(1,num_scaling_factors);
        plot_vec(1,:) = npads_mat(wire_res_ind,nind,:);
        plot_vec = plot_vec( plot_vec > 0);
        plot(1:length(plot_vec), plot_vec, 'color', colors{wire_res_ind}, 'linestyle', linestyles{nind})
    end
end
set(gca, 'xtick', 1:num_scaling_factors)
set(gca, 'xticklabel', node_labels)
colormap jet
xlabel('Process Node')
ylabel('Number of Power TSVs/Pads')
%set(gca,'yscale','log')
fixfigs(10,3,14,12)


