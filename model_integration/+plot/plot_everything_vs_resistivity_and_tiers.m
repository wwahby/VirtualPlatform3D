

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
    col_ind = mod(nind-1, length(colors))+1;
    plot(wire_resistivities*1e9, num_metal_levels_vec(nind,:),'color', colors{col_ind}, 'linestyle', '-');
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
    col_ind = mod(nind-1, length(colors))+1;
    plot(wire_resistivities*1e9, power_vec(nind,:), 'color', colors{col_ind}, 'linestyle', '-');
end
xlabel('Bulk Resistivity (\Omeganm)')
ylabel('Power (W)')
xlim([10 60])
fixfigs(2,3,14,12)

figure(3)
clf
hold on
for nind = 1:num_stacks
    col_ind = mod(nind-1, length(colors))+1;
    plot(wire_resistivities*1e9, npads_vec(nind,:), 'color', colors{col_ind}, 'linestyle', '-');
end
xlabel('Bulk Resistivity (\Omeganm)')
ylabel('Power Vias Per Tier')
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
    col_ind = mod(nind-1, length(colors))+1;
    area_per_die_mm2 = total_die_area_mm2/tiers(nind);
    total_tsv_area_allocation_mm2 = area_per_die_mm2 * tsv_max_area_fraction;
    
    plot(wire_resistivities*1e9, psn_tsv_area_vec_mm2(nind,:)/total_tsv_area_allocation_mm2, 'color', colors{col_ind}, 'linestyle', '-');
end
xlabel('Bulk Resistivity (\Omeganm)')
ylabel('Power Via Area / Max Via Area')
xlim([10 60])
set(gca,'yscale','log')
fixfigs(5,3,14,12)

figure(6)
clf
hold on
for nind = 1:num_stacks
    col_ind = mod(nind-1, length(colors))+1;
    area_per_die_mm2 = total_die_area_mm2/tiers(nind);
    total_tsv_area_allocation_mm2 = area_per_die_mm2;
    
    plot(wire_resistivities*1e9, psn_tsv_area_vec_mm2(nind,:)/total_tsv_area_allocation_mm2, 'color', colors{col_ind}, 'linestyle', '-');
end
xlabel('Bulk Resistivity (\Omeganm)')
ylabel('Power Via Area / Total Area')
xlim([10 60])
set(gca,'yscale','log')
fixfigs(6,3,14,12)


%% Cu
wire_res_ind = 1;
num_metal_levels_vec = zeros(num_stacks,num_scaling_factors);
power_vec = zeros(num_stacks,num_scaling_factors);
npads_vec = zeros(num_stacks,num_scaling_factors);
for nind = 1:num_stacks
    for scaling_ind = 1:num_scaling_factors
        num_metal_levels_vec(nind,scaling_ind) = num_metal_levels(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
        power_vec(nind, scaling_ind) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
        npads_vec(nind, scaling_ind) = npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);                                       
    end
end

num_metal_levels_vec( num_metal_levels_vec > 25) = NaN;
figure(7)
clf
hold on
surf(1:num_scaling_factors, 1:num_stacks, num_metal_levels_vec)
set(gca, 'xtick', 1:num_scaling_factors)
set(gca, 'xticklabel', node_labels)
colorbar
colormap jet
xlabel('Process Node')
ylabel('Number of Tiers')
caxis([2 11])
fixfigs(7,3,14,12)
%% W
wire_res_ind = 2;
num_metal_levels_vec = zeros(num_stacks,num_scaling_factors);
power_vec = zeros(num_stacks,num_scaling_factors);
npads_vec = zeros(num_stacks,num_scaling_factors);
for nind = 1:num_stacks
    for scaling_ind = 1:num_scaling_factors
        num_metal_levels_vec(nind,scaling_ind) = num_metal_levels(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
        power_vec(nind, scaling_ind) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
        npads_vec(nind, scaling_ind) = npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);                                       
    end
end

num_metal_levels_vec( num_metal_levels_vec > 25) = NaN;
figure(8)
clf
hold on
imagesc(1:num_scaling_factors, 1:num_stacks, num_metal_levels_vec)
set(gca, 'xtick', 1:num_scaling_factors)
set(gca, 'xticklabel', node_labels)
axis tight
colorbar
colormap jet
xlabel('Process Node')
ylabel('Number of Tiers')
caxis([2 11])
fixfigs(8,3,14,12)

%% Cu
wire_res_ind = 1;
num_metal_levels_vec = zeros(num_stacks,num_scaling_factors);
power_vec = zeros(num_stacks,num_scaling_factors);
npads_vec = zeros(num_stacks,num_scaling_factors);
for nind = 1:num_stacks
    for scaling_ind = 1:num_scaling_factors
        num_metal_levels_vec(nind,scaling_ind) = num_metal_levels(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
        power_vec(nind, scaling_ind) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
        npads_vec(nind, scaling_ind) = npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);                                       
    end
end

num_metal_levels_vec( num_metal_levels_vec > 25) = NaN;
num_metal_levels_vec = [num_metal_levels_vec NaN*ones(num_stacks,1) ; NaN*ones(1,num_scaling_factors + 1)];
figure(9)
clf
hold on
pcolor(1:num_scaling_factors+1, 1:num_stacks+1, num_metal_levels_vec)
set(gca, 'xtick', 1.5:1:num_scaling_factors+1)
set(gca, 'xticklabel', node_labels)
set(gca, 'ytick', 1.5:1:num_stacks+1)
set(gca, 'yticklabel', 1:num_stacks)
set(gca,'TickLength',[ 0 0 ])
axis tight
colorbar
colormap jet
xlabel('Process Node')
ylabel('Number of Tiers')
caxis([2 11])
fixfigs(9,3,14,12)

%% W
wire_res_ind = 2;
num_metal_levels_vec = zeros(num_stacks,num_scaling_factors);
power_vec = zeros(num_stacks,num_scaling_factors);
npads_vec = zeros(num_stacks,num_scaling_factors);
for nind = 1:num_stacks
    for scaling_ind = 1:num_scaling_factors
        num_metal_levels_vec(nind,scaling_ind) = num_metal_levels(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
        power_vec(nind, scaling_ind) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
        npads_vec(nind, scaling_ind) = npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);                                       
    end
end

num_metal_levels_vec( num_metal_levels_vec > 25) = NaN;
%num_metal_levels_vec = [num_metal_levels_vec NaN*ones(num_stacks,1) ; NaN*ones(1,num_scaling_factors + 1)];

figure(10)
clf
hold on
%pcolor(1:num_scaling_factors+1, 1:num_stacks+1, num_metal_levels_vec)
fig = surf_centered(1:num_scaling_factors, 1:num_stacks, num_metal_levels_vec);
set(gca, 'xtick', 1.5:1:num_scaling_factors+1)
set(gca, 'xticklabel', node_labels)
set(gca, 'ytick', 1.5:1:num_stacks+1)
set(gca, 'yticklabel', 1:num_stacks)
set(gca,'TickLength',[ 0 0 ])
axis tight
colorbar
colormap jet
xlabel('Process Node')
ylabel('Number of Tiers')
caxis([2 11])
fixfigs(10,3,14,12)

%% Transient power noise
% Cu
wire_res_ind = 2;
num_metal_levels_vec = zeros(num_stacks,num_scaling_factors);
power_vec = zeros(num_stacks,num_scaling_factors);
npads_vec = zeros(num_stacks,num_scaling_factors);
Npts_pow_transient = Npts_psn_tvec;
power_transient_mat = zeros(num_stacks, num_scaling_factors, Npts_pow_transient);
max_metal_levels_routable = 20;

tdecay_mat = zeros(num_stacks, num_scaling_factors);

figure(11)
clf
hold on
for nind = 1:num_stacks
    for scaling_ind = 1:num_scaling_factors
        num_metal_levels_vec(nind,scaling_ind) = num_metal_levels(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
        power_vec(nind, scaling_ind) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
        npads_vec(nind, scaling_ind) = npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
        if(routable_design(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind))
            power_transient_mat(nind, scaling_ind, :) = psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind}.output_cell{9};
            tvec = psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind}.output_cell{10};
            tstart = find(tvec  >= 100e-9, 1, 'first');
            tstop = find(tvec <= 130e-9, 1, 'last');

            tnorm_start = find(tvec >=500e-9, 1, 'first');
            tnorm_stop = find(tvec <= 800e-9, 1, 'last');
            final_rms_val = - rms( squeeze( power_transient_mat(nind, scaling_ind, tnorm_start:tnorm_stop) ) );
            [max_val, max_ind] = min( power_transient_mat(nind, scaling_ind, tstart:tstop) );
            max_delta = abs(max_val - final_rms_val);
            decay_threshold_ratio = 0.10; % percentage of max delta at which we'll say ringdown is complete
            decay_threshold = decay_threshold_ratio*max_delta;
            decay_diff_vec = abs( squeeze( power_transient_mat(nind, scaling_ind, :) ) - final_rms_val );
            if (abs((max_val - final_rms_val)/final_rms_val) < decay_threshold_ratio)
                tdecay_ind = tstart + max_ind;
            else
                tdecay_ind = find(decay_diff_vec(tstart:tstop) < decay_threshold, 1, 'last');
            end
            tdecay = tvec(tstart + tdecay_ind);
            tdecay_mat(nind, scaling_ind) = tdecay;
            
            
            plot(tvec(tstart:tstop)*1e9, squeeze(power_transient_mat(nind, scaling_ind, tstart:tstop))./rms(squeeze(power_transient_mat(nind, scaling_ind, tnorm_start:tnorm_stop))) )
        end

    end
end
num_metal_levels_vec( num_metal_levels_vec > max_metal_levels_routable) = NaN;

tdecay_mat(tdecay_mat == 0) = NaN;
figure(13)
clf
hold on
for nind = 2:num_stacks
    plot(tdecay_mat(nind, :)*1e9);
end
set(gca, 'xtick', 1:length(scaling_factors));
set(gca, 'xticklabels', node_labels);
xlabel('Node')
fixfigs(13,3,14,12)
ylabel('Transient decay time (ns)')



