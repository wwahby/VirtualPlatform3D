
colors = {'k', 'b', 'r', 'g'};
colors_dashed = {'k--', 'b--', 'r--', 'g--'};
colors_dot = {'k:', 'b:', 'r:', 'g:'};


tsv_num_vec = zeros(1,num_stacks);
twl_vec = zeros(1,num_stacks);

%fft
% panth_tsv_num = [0, 1050, 1921, 2475];
% panth_twl = [4927746, 4754600, 4745069, 4759862 ]; % um
% panth_twl_m = panth_twl/1e6;

%mult
panth_tsv_num = [0, 48513, 79682, 102994];
panth_twl = [29444308, 26169716, 23746536, 22470562 ]; % um
panth_twl_m = panth_twl/1e6;

for nind = 1:num_stacks
    tsv_num_vec(nind) = sum(tsv_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind,num_gates_ind}.per_layer);
    Iidf = ild_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind,num_gates_ind};
    gate_pitch_m = chip_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind,num_gates_ind}.gate_pitch;
    twl_vec(nind) = get_total_length(Iidf)*(gate_pitch_m);
end

twl_vec_m = twl_vec;

figure(1)
clf
hold on
grid on
plot(panth_tsv_num, 'k', 'linewidth', 1.5)
plot(tsv_num_vec, 'b', 'linewidth', 1.5)
set(gca, 'xtick', tiers)
xlabel('Number of Tiers')
ylabel('Number of MIVs')


figure(2)
clf
hold on
grid on
plot(panth_twl_m,'k', 'linewidth', 1.5)
plot(twl_vec_m, 'b', 'linewidth', 1.5)
set(gca, 'xtick', tiers)
xlabel('Number of Tiers')
ylabel('Total Wirelength (m)')