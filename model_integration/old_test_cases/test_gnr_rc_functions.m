clear all
close all

gnr_width = 7.5e-9; % (m)
gnr_length = 14e-6; % (m)
mfp_eff = 1e-6; % (m)
num_layers = 5; % (-)
rho_interlayer_cm = [0.3, 3, 30]; % (Ohm*cm)
Ef = 0.2; % (eV)
temp_K = 300; % (K)

contact_resistance = 0;
prob_backscattering = 0;
defect_mfp = mfp_eff;

rho_vec_m = rho_interlayer_cm/1e2;
max_segs = 300;

num_resistivities = length(rho_vec_m);
R_tc_mat = zeros(num_resistivities, max_segs);
R_k_mat = zeros(num_resistivities, max_segs);


for r_ind = 1:length(rho_vec_m)
    rho_interlayer = rho_vec_m(r_ind);
    R_tc_vec = zeros(1,max_segs);
    R_k_vec = zeros(1,max_segs);
    
    for N_segs = 1:max_segs
        [ R_tc, R_tc_pul, R_sc, R_sc_pul ] = xcm.gnr_calc_resistance_ns2013(gnr_width, gnr_length, mfp_eff, num_layers, rho_interlayer, Ef, temp_K, N_segs);
        [Rnet, Reff] = xcm.calc_gnr_resistance_kumar(num_layers, gnr_width, gnr_length, temp_K, defect_mfp, rho_interlayer, prob_backscattering, Ef, contact_resistance, N_segs);

        R_tc_vec(N_segs) = R_tc;
        R_k_vec(N_segs) = Reff;
    end
    
    R_tc_mat(r_ind,:) = R_tc_vec;
    R_k_mat(r_ind,:) = R_k_vec;
end
ymax = 1.2*max(max(R_tc_mat))/1e3;

%% Plots
colors = [0 0 1; 0 1 0; 1 0 0];
figure(1)
clf
hold on
for r_ind = 1:num_resistivities
    c_ind = mod(r_ind, 3) + 1;
    
    plot(R_k_mat(r_ind,:)/1e3, 'color', colors(c_ind,:), 'linestyle', '-')
    plot(R_tc_mat(r_ind,:)/1e3, 'color', colors(c_ind,:), 'linestyle', '--')
end
ylim([0, ymax])
xlabel('Number of segments')
ylabel('Resistance (k\Omega)')
fixfigs(1,2,14,12)


%
% figure(2)
% clf
% hold on
% for r_ind = 1:num_resistivities
%     plot(R_tc_mat(r_ind,:)/1e3)
% end
% xlabel('Number of segments')
% ylabel('Resistance (k\Omega)')
% fixfigs(2,2,14,12)
% 
% figure(3)
% clf
% hold on
% for r_ind = 1:num_resistivities
%     plot(R_k_mat(r_ind,:)/1e3)
% end
% ylim([0, ymax])
% xlabel('Number of segments')
% ylabel('Resistance (k\Omega)')
% fixfigs(3,2,14,12)


%% Resistance vs number of layers for a fixed number of segments

N_segs = 50;
max_layers = 10;
R_tc_layer_mat = zeros(num_resistivities, max_layers);
R_k_layer_mat = zeros(num_resistivities, max_layers);


for r_ind = 1:length(rho_vec_m)
    rho_interlayer = rho_vec_m(r_ind);
    R_tc_vec = zeros(1,max_layers);
    R_k_vec = zeros(1,max_layers);
    
    for num_layers = 1:max_layers
        [ R_tc, R_tc_pul, R_sc, R_sc_pul ] = xcm.gnr_calc_resistance_ns2013(gnr_width, gnr_length, mfp_eff, num_layers, rho_interlayer, Ef, temp_K, N_segs);
        [Rnet, Reff] = xcm.calc_gnr_resistance_kumar(num_layers, gnr_width, gnr_length, temp_K, defect_mfp, rho_interlayer, prob_backscattering, Ef, contact_resistance, N_segs);

        R_tc_vec(num_layers) = R_tc;
        R_k_vec(num_layers) = Reff;
    end
    
    R_tc_layer_mat(r_ind,:) = R_tc_vec;
    R_k_layer_mat(r_ind,:) = R_k_vec;
end
ymax = 1.2*max(max(R_tc_layer_mat))/1e3;

colors = [0 0 1; 0 1 0; 1 0 0];
figure(2)
clf
hold on
for r_ind = 1:num_resistivities
    c_ind = mod(r_ind-1, 3)+1;
    
    plot(R_k_layer_mat(r_ind,:)/1e3, 'color', colors(c_ind,:), 'linestyle', '-')
    plot(R_tc_layer_mat(r_ind,:)/1e3, 'color', colors(c_ind,:), 'linestyle', '--')
end
ylim([0, ymax])
xlabel('Number of layers')
ylabel('Resistance (k\Omega)')
fixfigs(2,2,14,12)


%% Reproduce Vachan's figure

N_segs = 50; % NS2013 formula is using N_segs = 4.
max_layers = 10;
gnr_width = 7.5e-9;
gnr_length_vec = [7e-6, 14e-6];
rho_interlayer_cm = 30;
rho_interlayer = rho_interlayer_cm/1e2;

num_lengths = length(gnr_length_vec);
R_tc_layer_mat = zeros(num_lengths, max_layers);
R_k_layer_mat = zeros(num_lengths, max_layers);

for l_ind = 1:num_lengths
    
    gnr_length = gnr_length_vec(l_ind);
    
    R_tc_vec = zeros(1,max_layers);
    R_k_vec = zeros(1,max_layers);
    
    for num_layers = 1:max_layers
        % Note that ns2013 is most accurate for N_seg = 4.
        [ R_tc, R_tc_pul, R_sc, R_sc_pul ] = xcm.gnr_calc_resistance_ns2013(gnr_width, gnr_length, mfp_eff, num_layers, rho_interlayer, Ef, temp_K, 4);
        [Rnet, Reff] = xcm.calc_gnr_resistance_kumar(num_layers, gnr_width, gnr_length, temp_K, defect_mfp, rho_interlayer, prob_backscattering, Ef, contact_resistance, N_segs);

        R_tc_vec(num_layers) = R_tc;
        R_k_vec(num_layers) = Reff;
    end
    
    R_tc_layer_mat(l_ind,:) = R_tc_vec;
    R_k_layer_mat(l_ind,:) = R_k_vec;
end

ymax = 1.2*max(max(R_tc_layer_mat))/1e3;

colors = [0 0 1; 0 1 0; 1 0 0];
figure(3)
clf
hold on
for l_ind = 1:num_lengths
    c_ind = mod(l_ind-1, 3)+1;
    
    plot(R_k_layer_mat(l_ind,:)/1e3, 'color', colors(c_ind,:), 'linestyle', '-')
    plot(R_tc_layer_mat(l_ind,:)/1e3, 'color', colors(c_ind,:), 'linestyle', '--')
end
ylim([0, ymax])
xlabel('Number of layers')
ylabel('Resistance (k\Omega)')
fixfigs(3,2,14,12)
