%function plot_sweep_data( sweep, sweep_data, simulation) % just use as a
%script for now

%% Unpack sweep data

num_stacks = length(tiers);
num_perms = length(rel_permittivities);
num_thicks = length(thicknesses);

if (simulation.freq_binsearch == 1)
    num_freqs = 1; % skip freq loops
else
    num_freqs = length(frequencies);
end
num_cooling_configs = length(heat_fluxes);
num_decaps = length(decap_ratios);
num_wire_resistivities = length(wire_resistivities);
num_wire_flags = length(wire_material_flags);
num_scaling_factors = length(scaling_factors);
num_barrier_thicknesses = length(barrier_thicknesses);
num_barrier_resistivities = length(barrier_resistivities);

cind = num_cooling_configs;
dind = num_decaps;
thind = num_thicks;
nind = num_stacks;
pind = num_perms;
freq_ind = num_freqs;
wire_res_ind = num_wire_resistivities;
wire_flag_ind = num_wire_flags;
scaling_ind = num_scaling_factors;
bar_thick_ind = num_barrier_thicknesses;
bar_res_ind = num_barrier_resistivities;

%% Plot options
color_vec = [0 0 0; 0 0 1; 0 1 0; 1 0 0];
linewidth = 2;
label_size = 14;
axis_font_size = 12;

%% Plot WLD (GP) for each configuration
fignum = 1;
figure(fignum)
clf
hold on
max_wires = 0;
for nind = 1:num_stacks
    WLD = ild_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind};
    max_i = max(WLD);
    max_wires = max(max_wires, max_i);
    
    color_ind = mod(nind-1, 4)+1;
    plot(WLD, 'color', color_vec(color_ind, :), 'linewidth', linewidth)
end
set(gca,'yscale','log')
set(gca,'xscale','log')
ylim([1e-1, max_wires*10])
xlabel('Wire Length (GP)')
ylabel('Number of wires')
fixfigs(fignum, linewidth, label_size, axis_font_size)


%% Plot WLD (um) for each configuration
fignum = fignum+1;
figure(fignum)
clf
hold on
max_wires = 0;
for nind = 1:num_stacks
    WLD = ild_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind};
    max_i = max(WLD);
    max_wires = max(max_wires, max_i);
    
    color_ind = mod(nind-1, 4)+1;
    xvec = (1:length(WLD)) * design.gate_pitch*1e6;
    plot(xvec, WLD, 'color', color_vec(color_ind, :), 'linewidth', linewidth)
end
set(gca,'yscale','log')
set(gca,'xscale','log')
ylim([1e-1, max_wires*10])
xlabel('Wire Length (um)')
ylabel('Number of wires')
fixfigs(fignum, linewidth, label_size, axis_font_size)

%% Plot Total WL (um) for each configuration

mult_256_skl_miv_wl = [30224686, 26169716, 23746536, 22470532];
cf_fft_skl_miv_wl = [4927746, 4754600, 4745069, 4759862];
cf_rca_skl_miv_wl = [1578160, 1499774, 1466258, 1451201];
des_perf_skl_miv_wl = [566977, 478166, 432728, 415356];

fignum = fignum+1;
figure(fignum)
clf
hold on
total_wl = zeros(1,num_stacks);
for nind = 1:num_stacks
    color_ind = mod(nind-1, 4)+1;
    
    WLD = ild_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind};
    xvec = (1:length(WLD)) * design.gate_pitch*1e6;
    
    total_wl(nind) = sum(WLD.*xvec);
end

bar(tiers, [total_wl; des_perf_skl_miv_wl]', 1)
set(gca,'xtick', tiers)
xlabel('Number of tiers')
ylabel('Total WL (m)')
fixfigs(fignum, linewidth, label_size, axis_font_size)
    
                                          




