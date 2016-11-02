%function plot_sweep_data( sweep, sweep_data, simulation) % just use as a
%script for now
fignum = 1;

%% EPC, Freq, Temp, Total Power, Total Power Density

colors = [0 0 0; 0 0 1; 0 1 0; 1 0 0];
freq_mat = zeros(num_stacks, num_cooling_configs*num_thicks);
power_mat = zeros(num_stacks, num_cooling_configs*num_thicks);
power_density_mat = zeros(num_stacks, num_cooling_configs*num_thicks);
wire_power_mat = zeros(num_stacks, num_cooling_configs*num_thicks);
rep_power_mat = zeros(num_stacks, num_cooling_configs*num_thicks);
dynamic_power_mat = zeros(num_stacks, num_cooling_configs*num_thicks);
leakage_power_mat = zeros(num_stacks, num_cooling_configs*num_thicks);
area_mat = zeros(num_stacks, num_cooling_configs*num_thicks);
temperature_mat = zeros(num_stacks, num_cooling_configs*num_thicks);

A = zeros(num_stacks, num_cooling_configs*num_thicks);
B = zeros(num_stacks, num_cooling_configs*num_thicks);
C = zeros(num_stacks, num_cooling_configs*num_thicks);

for thind = 1:num_thicks
    for nind = 1:num_stacks
        for cind = 1:num_cooling_configs
            A(nind, (cind-1)*num_thicks + thind ) = cind;
            B(nind, (cind-1)*num_thicks + thind ) = nind;
            C(nind, (cind-1)*num_thicks + thind ) = thind;
            
            freq_mat(nind, (cind-1)*num_thicks + thind ) = freq(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
            power_mat(nind, (cind-1)*num_thicks + thind ) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
            power_density_mat(nind, (cind-1)*num_thicks + thind ) = power_density(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
            wire_power_mat(nind, (cind-1)*num_thicks + thind ) = wire_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
            rep_power_mat(nind, (cind-1)*num_thicks + thind ) = rep_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
            dynamic_power_mat(nind, (cind-1)*num_thicks + thind ) = dynamic_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
            leakage_power_mat(nind, (cind-1)*num_thicks + thind ) = leakage_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
            area_mat(nind, (cind-1)*num_thicks + thind ) = chip_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind}.area_per_layer_m2;                               
            temperature_mat(nind, (cind-1)*num_thicks + thind ) = temperature(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
        end
    end
end
Tclk_mat = 1./freq_mat;
epc_mat = power_mat .* Tclk_mat;
epc_mat = epc_mat*1e9;

on_chip_power_frac_mat = (wire_power_mat + rep_power_mat)./power_mat;



%%
tier_cell = {};
for nind = 1:num_stacks
    tier_cell{nind} =  num2str(tiers(nind));
end
colors = {'b', 'y', 'g', 'r'};

figure(1)
clf
b = bar(freq_mat/1e9, 1, 'grouped');
% for bbb = 1:length(b)
%     b(bbb).FaceColor = colors{bbb};
% end
colormap(jet)
grid on
xlabel('Tiers')
set(gca,'xticklabel',tier_cell)
b(1).FaceColor = 'b';
b(2).FaceColor ='y';
b(3).FaceColor = 'r';
ylabel('Maximum Frequency (GHz)')
fixfigs(1,3,14,12)

figure(11)
clf
b = bar(freq_mat/freq_mat(1,1), 1, 'grouped');
% for bbb = 1:length(b)
%     b(bbb).FaceColor = colors{bbb};
% end
colormap(jet)
grid on
xlabel('Tiers')
set(gca,'xticklabel',tier_cell)
ylabel('Relative maximum clock rate')
b(1).FaceColor = 'b';
b(2).FaceColor ='y';
b(3).FaceColor = 'r';
fixfigs(11,3,14,12)
%%
figure(2)
clf
b = bar(epc_mat, 1, 'grouped');
% for bbb = 1:length(b)
%     b(bbb).FaceColor = colors{bbb};
% end
colormap(jet)
xlabel('Tiers')
set(gca,'xticklabel',tier_cell)
ylabel('Energy per Cycle (nJ)')
fixfigs(2,3,14,12)

figure(3)
clf
b = bar(power_mat, 1, 'grouped');
% for bbb = 1:length(b)
%     b(bbb).FaceColor = colors{bbb};
% end
colormap(jet)
xlabel('Tiers')
set(gca,'xticklabel',tier_cell)
ylabel('Power (W)')
fixfigs(3,3,14,12)

figure(4)
clf
b = bar(temperature_mat, 1, 'grouped');
% for bbb = 1:length(b)
%     b(bbb).FaceColor = colors{bbb};
% end
colormap(jet)
%ylim([0 90])
xlabel('Tiers')
set(gca,'xticklabel',tier_cell)
ylabel('Temperature (^\circC)')
fixfigs(4,3,14,12)

colors = {'b', 'r', 'g', 'y'};
linestyles = {'-', '--'};
figure(5)
clf
hold on
for thind = 1:num_thicks
    linestyle = linestyles{thind};
    for cind = 1:num_cooling_configs
        plot(tiers, freq_mat(:,(cind-1)*num_thicks + thind)/1e9, 'color', colors{cind}, 'linestyle', linestyle)
    end
end
xlabel('Tiers')
grid on
set(gca,'xtick',1:num_stacks)
ylabel('Maximum Frequency (GHz)')
fixfigs(5,3,14,12)

colors = {'b', 'r', 'g', 'y', 'm', 'c', 'k'};
figure(6)
clf
hold on
for thind = 1:num_thicks
    linestyle = linestyles{thind};
    for cind = 1:num_cooling_configs
        plot(tiers, on_chip_power_frac_mat(:,(cind-1)*num_thicks + thind),  'color', colors{cind}, 'linestyle', linestyle)
    end
end
xlabel('Tiers')
ylim([0.5 0.85])
grid on
set(gca,'xtick',1:num_stacks)
ylabel('On-Chip Communication Power Fraction')
fixfigs(6,3,14,12)

colors = {'b', 'r', 'g', 'y'};
figure(7)
clf
hold on
for thind = 1:num_thicks
    linestyle = linestyles{thind};
    for cind = 1:num_cooling_configs
        plot(tiers, leakage_power_mat(:,(cind-1)*num_thicks + thind),  'color', colors{cind}, 'linestyle', linestyle)
    end
end
xlabel('Tiers')
%ylim([0.5 0.85])
grid on
set(gca,'xtick',1:num_stacks)
ylabel('Leakage Power (W)')
fixfigs(7,3,14,12)

colors = {'b', 'r', 'g', 'y'};
figure(8)
clf
hold on
for thind = 1:num_thicks
    linestyle = linestyles{thind};
    for cind = 1:num_cooling_configs
        plot(tiers, epc_mat(:,(cind-1)*num_thicks + thind),  'color', colors{cind}, 'linestyle', linestyle)
    end
end
xlabel('Tiers')
%ylim([0.5 0.85])
grid on
set(gca,'xtick',1:num_stacks)
ylabel('Energy Per Cycle (nJ)')
fixfigs(8,3,14,12)