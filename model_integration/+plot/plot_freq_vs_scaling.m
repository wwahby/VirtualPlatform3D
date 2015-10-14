%function plot_sweep_data( sweep, sweep_data, simulation) % just use as a
%script for now
fignum = 1;

%% EPC vs tiers

colors = [0 0 0; 0 0 1; 0 1 0; 1 0 0];
freq_mat = zeros(num_stacks, num_scaling_factors);
power_mat = zeros(num_stacks, num_scaling_factors);
power_density_mat = zeros(num_stacks, num_scaling_factors);
wire_power_mat = zeros(num_stacks, num_scaling_factors);
rep_power_mat = zeros(num_stacks, num_scaling_factors);
dynamic_power_mat = zeros(num_stacks, num_scaling_factors);
leakage_power_mat = zeros(num_stacks, num_scaling_factors);
area_mat = zeros(num_stacks, num_scaling_factors);

for nind = 1:num_stacks
    for scaling_ind = 1:num_scaling_factors
        freq_mat(nind, scaling_ind) = freq(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
        power_mat(nind, scaling_ind) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
        power_density_mat(nind, scaling_ind) = power_density(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
        wire_power_mat(nind, scaling_ind) = wire_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
        rep_power_mat(nind, scaling_ind) = rep_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
        dynamic_power_mat(nind, scaling_ind) = dynamic_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
        leakage_power_mat(nind, scaling_ind) = leakage_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
        area_mat(nind, scaling_ind) = chip_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind}.area_per_layer_m2;                               
    end
end
Tclk_mat = 1./freq_mat;
epc_mat = power_mat .* Tclk_mat;
epc_mat = epc_mat'*1e9;

figure(1)
clf
hold on
for nind = 1:num_stacks
    plot(1:num_scaling_factors, freq_mat(nind,:)/1e9, 'color', colors(nind,:))
end
set(gca,'Xtick',1:num_scaling_factors)
set(gca,'XtickLabel', node_labels )
xlabel('Process Node')
ylabel('Maximum Frequency (GHz)')
fixfigs(1,3,14,12)

figure(2)
clf
hold on
for nind = 1:num_stacks
    plot(1:num_scaling_factors, epc_mat(:,nind), 'color', colors(nind,:))
end
set(gca,'Xtick',1:num_scaling_factors)
set(gca,'XtickLabel', node_labels )
xlabel('Process Node')
ylabel('Energy Per Cycle (nJ)')
fixfigs(2,3,14,12)

figure(3)
clf
hold on
for nind = 1:num_stacks
    plot(1:num_scaling_factors, power_mat(nind,:), 'color', colors(nind,:))
end
set(gca,'Xtick',1:num_scaling_factors)
set(gca,'XtickLabel', node_labels )
xlabel('Process Node')
ylabel('Power (W)')
fixfigs(3,3,14,12)

figure(4)
clf
hold on
for nind = 1:num_stacks
    plot(1:num_scaling_factors, power_mat(nind,:), 'color', colors(nind,:))
end
set(gca,'Xtick',1:num_scaling_factors)
set(gca,'XtickLabel', node_labels )
xlabel('Process Node')
ylabel('Power (W)')
fixfigs(4,3,14,12)

figure(5)
clf
hold on
for nind = 1:num_stacks
    plot(1:num_scaling_factors, power_density_mat(nind,:)/1e4, 'color', colors(nind,:))
end
set(gca,'Xtick',1:num_scaling_factors)
set(gca,'XtickLabel', node_labels )
xlabel('Process Node')
ylabel('Power Density (W/cm^2)')
fixfigs(5,3,14,12)

fignum = 5;
fignum = fignum + 1;
figure(fignum)
clf
hold on
for nind = 1:num_stacks
    plot(1:num_scaling_factors, wire_power_mat(nind,:), 'color', colors(nind,:))
end
set(gca,'Xtick',1:num_scaling_factors)
set(gca,'XtickLabel', node_labels )
xlabel('Process Node')
ylabel('Wire Power (W)')

fignum = fignum + 1;
figure(fignum)
clf
hold on
for nind = 1:num_stacks
    plot(1:num_scaling_factors, rep_power_mat(nind,:), 'color', colors(nind,:))
end
set(gca,'Xtick',1:num_scaling_factors)
set(gca,'XtickLabel', node_labels )
xlabel('Process Node')
ylabel('Repeater Power (W)')

fignum = fignum + 1;
figure(fignum)
clf
hold on
for nind = 1:num_stacks
    plot(1:num_scaling_factors, dynamic_power_mat(nind,:), 'color', colors(nind,:))
end
set(gca,'Xtick',1:num_scaling_factors)
set(gca,'XtickLabel', node_labels )
xlabel('Process Node')
ylabel('Dynamic Power (W)')

fignum = fignum + 1;
figure(fignum)
clf
hold on
for nind = 1:num_stacks
    plot(1:num_scaling_factors, leakage_power_mat(nind,:), 'color', colors(nind,:))
end
set(gca,'Xtick',1:num_scaling_factors)
set(gca,'XtickLabel', node_labels )
xlabel('Process Node')
ylabel('Leakage Power (W)')

fixfigs(6:fignum, 3, 14, 12)

%%
area_mat_cm2 = area_mat*1e4;
fignum_pre = fignum + 1;
fignum = fignum + 1;
figure(fignum)
clf
hold on
for nind = 1:num_stacks
    plot(1:num_scaling_factors, wire_power_mat(nind,:)./area_mat_cm2(nind,:), 'color', colors(nind,:))
end
set(gca,'Xtick',1:num_scaling_factors)
set(gca,'XtickLabel', node_labels )
xlabel('Process Node')
ylabel('Wire Power Density (W/cm^2)')

fignum = fignum + 1;
figure(fignum)
clf
hold on
for nind = 1:num_stacks
    plot(1:num_scaling_factors, rep_power_mat(nind,:)./area_mat_cm2(nind,:), 'color', colors(nind,:))
end
set(gca,'Xtick',1:num_scaling_factors)
set(gca,'XtickLabel', node_labels )
xlabel('Process Node')
ylabel('Repeater Power Density (W/cm^2)')

fignum = fignum + 1;
figure(fignum)
clf
hold on
for nind = 1:num_stacks
    plot(1:num_scaling_factors, dynamic_power_mat(nind,:)./area_mat_cm2(nind,:), 'color', colors(nind,:))
end
set(gca,'Xtick',1:num_scaling_factors)
set(gca,'XtickLabel', node_labels )
xlabel('Process Node')
ylabel('Dynamic Power Density (W/cm^2)')

fignum = fignum + 1;
figure(fignum)
clf
hold on
for nind = 1:num_stacks
    plot(1:num_scaling_factors, leakage_power_mat(nind,:)./area_mat_cm2(nind,:), 'color', colors(nind,:))
end
set(gca,'Xtick',1:num_scaling_factors)
set(gca,'XtickLabel', node_labels )
xlabel('Process Node')
ylabel('Leakage Power Density (W/cm^2)')

fixfigs(fignum_pre:fignum, 3, 14, 12)




