%function plot_sweep_data( sweep, sweep_data, simulation) % just use as a
%script for now

%% Temp

figure(1)
clf
hold on
colors = [0 0 0; 0 0 1; 0 1 0; 1 0 0];
linestyles = {'-','--'};
for nind = 1:num_stacks
    tvec = zeros(1, num_forced_powers);
    for cind = 1:num_cooling_configs
        for forced_power_ind = 1:num_forced_powers
            tvec(forced_power_ind) = temperature(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind);
        end
        plot(power_forced_vec, tvec, 'color', colors(nind,:), 'linestyle', linestyles{cind} )
    end
    
end
ylim([20 100])
xlabel('Power (W)')
ylabel('Maximum Temperature (C)')
fixfigs(1,3,14,12)


