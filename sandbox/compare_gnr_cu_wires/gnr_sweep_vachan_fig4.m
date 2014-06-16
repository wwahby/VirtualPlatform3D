% Sweep mGNR params
clear all
close all

layer_vec = 1:1:10;
length_vec = (1:1:20)*1e-6;
num_lengths = length(length_vec);


temp = 300;
gnr_width = 7.5e-9;
rho_interlayer = 3e-1;
prob_backscattering = 0.0;
Ef = 0.2;
contact_resistance = 0;
defect_mfp = 1000e-9;

total_layers = length(layer_vec);
Rnet_mat = zeros(total_layers,num_lengths);
Reff_mat = zeros(total_layers,num_lengths);
R_nishad_top_mat = zeros(total_layers,num_lengths);
R_nishad_side_mat = zeros(total_layers,num_lengths);

height_dielectric = 2e-6;
epsrd = 16.5;

for nind = 1:length(layer_vec)
    for lind = 1:num_lengths
        
        num_layers = layer_vec(nind);
        gnr_length = length_vec(lind);

        [Rnet Reff] = calc_gnr_params(num_layers,gnr_width,gnr_length,temp,defect_mfp,rho_interlayer,prob_backscattering,Ef,contact_resistance);
        Rnet_mat(nind,lind) = Rnet;
        Reff_mat(nind,lind) = Reff;
        
        [delay_top delay_side R_top R_side L C] = calc_gnr_params_nishad(num_layers,gnr_width,gnr_length,temp,defect_mfp,rho_interlayer,prob_backscattering,Ef,contact_resistance,epsrd,height_dielectric);
        %[delay_top delay_side R_top R_side L C] = calc_gnr_params_combined(num_layers,gnr_width,gnr_length,temp,defect_mfp,rho_interlayer,prob_backscattering,Ef,contact_resistance,epsrd,height_dielectric);
        
        R_nishad_top_mat(nind,lind) = R_top;
        R_nishad_side_mat(nind,lind) = R_side;
    end
end


%%
figure(1)
clf
plot(length_vec,Rnet_mat)

figure(2)
clf
plot(length_vec,Reff_mat)



%% Compare to Vachan's paper
gnr_kumar_length = [1 5 10 15 20]*1e-6; % (m) GNR length
gnr_kumar_width = 10e-9; % (m) GNR width

gnr_kumar_layers = [2 4 6 8];
gnr_kumar_R = [13.66 20.05 27.53 34.86 42.26;
               13.66 19.74 24.73 28.62 32.52;
               13.66 19.74 24.18 26.83 29.40;
               13.66 19.74 24.10 26.29 28.16]*1e3;
           
gnr_my_kumar_R = [21.12 103.6 207.1 310.5 414.0]*1e3;
           
gnr_nishad_R = [8.11 40.41 80.81 121.2 161.6;
               4.13 20.20 40.41 61.61 80.82;
               2.85 13.50 26.95 40.42 53.88;
               2.24 10.15 20.23 30.32 40.42]*1e3;


figure(3)
clf
hold all
set(gca,'ColorOrder',[0 0 0; 0 0 1; 0 1 0; 1 0 0]);
set(gca,'LineStyleOrder','-')
% plot(gnr_kumar_length*1e6,gnr_kumar_R/1e3)
for pind = 1:length(gnr_kumar_layers)
    plot(gnr_kumar_length*1e6,gnr_kumar_R(pind,:)/1e3)
end

set(gca,'ColorOrder',[0 0 0; 0 0 1; 0 1 0; 1 0 0]);
set(gca,'LineStyleOrder','--')
% plot(length_vec*1e6,Reff_mat/1e3)
for pind = 1:length(layer_vec)
    plot(length_vec*1e6,Reff_mat(pind,:)/1e3)
end
fixfigs(3,3,14,12)

%% Compare Paper, Kumar model, and Nishad model

figure(4)
clf
hold all
% set(gca,'ColorOrder',[0 0 0; 0 0 1; 0 1 0; 1 0 0]);
% set(gca,'LineStyleOrder','-')
% % plot(gnr_kumar_length*1e6,gnr_kumar_R/1e3)
% for pind = 1:length(layer_vec)
%     plot(gnr_kumar_length*1e6,gnr_kumar_R(pind,:)/1e3)
% end

set(gca,'ColorOrder',[0 0 0; 0 0 1; 0 1 0; 1 0 0]);
set(gca,'LineStyleOrder','-')
% plot(length_vec*1e6,Reff_mat/1e3)
for pind = 1:length(layer_vec)
    plot(length_vec*1e6,Reff_mat(pind,:)/1e3)
end

set(gca,'ColorOrder',[0 0 0; 0 0 0.5; 0 0.5 0; 0.5 0 0]);
set(gca,'LineStyleOrder','-')
% plot(length_vec*1e6,Reff_mat/1e3)
for pind = 1:length(layer_vec)
    plot(length_vec*1e6,R_nishad_top_mat(pind,:)/1e3)
end
fixfigs(4,3,14,12)

%% Actual fig 4
gnr_kumar_num_layers = 1:10;
gnr_kumar_R_fig4_7um = [ 48.30508475	29.66101695	25.42372881	24.15254237	23.72881356	23.72881356	23.72881356	23.72881356	23.72881356	23.72881356 ];
gnr_kumar_R_fig4_14um = [ 95.33898305	52.54237288	39.83050847	34.3220339	31.77966102	30.08474576	29.66101695	29.66101695	29.23728814	28.81355932 ];
lind = 7;
R_kumar_layers_7um = Reff_mat(:,lind);
R_nishad_layers_7um = R_nishad_top_mat(:,lind);

lind = 14;
R_kumar_layers_14um = Reff_mat(:,lind);
R_nishad_layers_14um = R_nishad_top_mat(:,lind);

figure(5)
clf
plot(gnr_kumar_num_layers,gnr_kumar_R_fig4_7um,'k-')
hold on
plot(gnr_kumar_num_layers,gnr_kumar_R_fig4_14um,'k--')

plot(layer_vec,R_kumar_layers_7um/1e3,'b-')
plot(layer_vec,R_nishad_layers_7um/1e3,'r-')

plot(layer_vec,R_kumar_layers_14um/1e3,'b--')
plot(layer_vec,R_nishad_layers_14um/1e3,'r--')
ylim([0 100])
xlabel('Number of layers')
ylabel('Resistance (k\Omega)')
fixfigs(5,3,14,12)







