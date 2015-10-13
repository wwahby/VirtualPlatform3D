%function plot_sweep_data( sweep, sweep_data, simulation) % just use as a
%script for now
fignum = 1;
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


%% Heat transfer coefficient determined by sweep

h_vec = zeros(1,num_stacks);
for nind = 1:num_stacks
    h_vec(nind) = h_coeff(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind);
end

figure(1)
clf
plot(tiers, h_vec/1e4)
xlabel('Number of tiers')
ylabel('Heat Transfer Coefficient (W/cm^2K)')
fixfigs(1,3,14,12)

%% Actual tier power

pvec = zeros(1,num_stacks);
freq_vec = zeros(1, num_stacks);
for nind = 1:num_stacks
    pvec(nind) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind);
    freq_vec(nind) = freq(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind);
end

figure(2)
clf
plot(tiers, pvec)
xlabel('Number of tiers')
ylabel('Power (W)')
fixfigs(2,3,14,12)


