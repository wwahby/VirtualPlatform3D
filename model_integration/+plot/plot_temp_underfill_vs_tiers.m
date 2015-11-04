%function plot_sweep_data( sweep, sweep_data, simulation) % just use as a
%script for now
fignum = 1;

%% EPC, Freq, Temp, Total Power, Total Power Density

colors = {'k', 'b', 'g', 'r'};
linestyles = {'-', '--', ':', '-.'};




temperature_mat = zeros(num_stacks, num_thermal_conductivities);
for nind = 1:num_stacks
    for k_ind = 1:num_thermal_conductivities
        for scaling_ind = 1:num_scaling_factors
            if (routable_design(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind) == 1)
                temperature_mat(nind, k_ind) = temperature(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
            else
                temperature_mat(nind, k_ind) = 0;
            end  
        end
    end
end


figure(1)
clf
hold on
for nind = 1:1
    plot(thermal_conductivities, temperature_mat(nind,:), 'color', colors{nind})
end
xlabel('Underfill Thermal Conductivity (W/mK)')
ylabel('Maximum temperature (^\circC)')
fixfigs(1,3,14,12)