
num_metal_levels_vec_tsv = zeros(1,num_stacks);
power_vec_tsv = zeros(1,num_stacks);
npads_vec_tsv = zeros(1,num_stacks);
wire_res_ind = 1;
thind = 2;
for nind = 1:num_stacks
    num_metal_levels_vec_tsv(nind,wire_res_ind) = num_metal_levels(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
    power_vec_tsv(nind) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
    npads_vec_tsv(nind) = npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
end

colors = {'k', 'b', 'g', 'r'};

num_metal_levels_vec_m3d_w = zeros(1,num_stacks);
power_vec_m3d_w = zeros(1,num_stacks);
npads_vec_m3d_w = zeros(1,num_stacks);
wire_res_ind = 2;
thind = 1;
for nind = 1:num_stacks
    num_metal_levels_vec_m3d_w(nind,wire_res_ind) = num_metal_levels(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
    power_vec_m3d_w(nind) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
    npads_vec_m3d_w(nind) = npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
end

num_metal_levels_vec_m3d_cu = zeros(1,num_stacks);
power_vec_m3d_cu = zeros(1,num_stacks);
npads_vec_m3d_cu = zeros(1,num_stacks);
wire_res_ind = 1;
thind = 1;
for nind = 1:num_stacks
    num_metal_levels_vec_m3d_cu(nind,wire_res_ind) = num_metal_levels(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
    power_vec_m3d_cu(nind) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
    npads_vec_m3d_cu(nind) = npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind);
end

pow_bar = [ power_vec_tsv ; power_vec_m3d_w ; power_vec_m3d_cu ]';

pow_bar_norm = bsxfun(@rdivide, pow_bar(2:end, :), pow_bar(1,:) );
pow_bar_norm_m3dw = pow_bar_norm(:, 1:2);



% TSV vs M3DW and M3DC - Raw Power
figure(1)
clf
hold on
grid on
b = bar(pow_bar, 1, 'grouped');
set(gca,'xtick',1:4)
b(1).FaceColor = 'b';
b(2).FaceColor = 'y';
b(3).FaceColor = 'r';
xlabel('Number of Tiers')
ylabel('Power Consumption (W)')
fixfigs(1,2,14,12)

%TSV vs M3DW and M3DC - Norm power
figure(2)
clf
hold on
grid on
b = bar(pow_bar_norm, 1, 'grouped');
set(gca, 'xtick', 1:3)
set(gca, 'xticklabel', 2:4)
b(1).FaceColor = 'b';
b(2).FaceColor = 'y';
b(3).FaceColor = 'r';
xlabel('Number of Tiers')
ylabel('Power Consumption (Normalized to 2D case)')
fixfigs(2,2,14,12)

%TSV vs M3DW only - Norm power
figure(3)
clf
hold on
grid on
b = bar(pow_bar_norm_m3dw, 1, 'grouped');
set(gca, 'xtick', 1:3)
set(gca, 'xticklabel', 2:4)
b(1).FaceColor = 'b';
b(2).FaceColor = 'y';
% b(3).FaceColor = 'r';
xlabel('Number of Tiers')
ylabel('Power Consumption (Normalized to 2D)')
fixfigs(3,2,14,12)