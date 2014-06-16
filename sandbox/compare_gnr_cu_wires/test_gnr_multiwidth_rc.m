%% Test the new combined R/C/delay calculation
clear all
close all

layer_vec = 5;
length_vec = (100)*1e-6;
num_lengths = length(length_vec);


rho_interlayer = 3e-1;
prob_backscattering = 0.0;
Ef = 0.2;
contact_resistance = 0;
mfp_defect = 1000e-9;

gnr_widths = (1:100)*1e-9;
temp_K = 300;
epsrd = 3;
height_dielectric = 300e-9;


for nind = 1:length(layer_vec)
    for lind = 1:num_lengths
        
        num_layers = layer_vec(nind);
        gnr_length = length_vec(lind);
        
        [delay_top_vec delay_side_vec R_top_vec R_top_alt_vec R_side_vec L_vec C_gnr_vec C_gnr_raw_vec Nch_vec mfp_eff_vec] = ...
            calc_gnr_params_combined_multiple_widths( ...
            num_layers, gnr_widths, gnr_length, temp_K, mfp_defect, ...
            rho_interlayer, prob_backscattering, Ef,contact_resistance, epsrd, height_dielectric );

    end
end

for wind = 1:length(gnr_widths)
    gnr_width = gnr_widths(wind);
    [delay_top delay_side R_top R_side L C] = calc_gnr_params_combined( ...
                num_layers, gnr_width, gnr_length, temp_K, mfp_defect, rho_interlayer, prob_backscattering, ...
                Ef, contact_resistance, epsrd, height_dielectric );
    R_top_k(wind) = R_top;
    R_side_k(wind) = R_side;
end
%% Delay

pl_widths = gnr_widths*1e9;

figure(1)
clf
semilogy(pl_widths, delay_top_vec*1e9, 'b')
hold on
semilogy(pl_widths, delay_side_vec*1e9, 'r')
xlim([min(pl_widths) max(pl_widths)])
xlabel('GNR width (nm)')
ylabel('Delay (ns)')
fixfigs(1,3,14,12)

%% Resistance
figure(2)
clf
plot(pl_widths, R_top_vec/1e3, 'b')
hold on
plot(pl_widths, R_top_alt_vec/1e3, 'c')
plot(pl_widths, R_side_vec/1e3, 'r')
plot(pl_widths, R_top_k/1e3, 'b--')
plot(pl_widths, R_side_k/1e3, 'r--')
set(gca,'yscale','log')
xlim([min(pl_widths) max(pl_widths)])
xlabel('GNR width (nm)')
ylabel('Resistance (k\Omega)')
fixfigs(2,3,14,12)

%% Capacitance
figure(3)
clf
plot(pl_widths, C_gnr_vec*1e15, 'b')
hold on
plot(pl_widths, C_gnr_raw_vec*1e15, 'r')
set(gca,'yscale','log')
xlim([min(pl_widths) max(pl_widths)])
ylim([1e-1 1e1])
xlabel('GNR width (nm)')
ylabel('Capacitance (fF)')
fixfigs(3,3,14,12)