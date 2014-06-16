% Sweep mGNR params

layer_vec = 1:1:10;
length_vec = (7)*1e-6;
num_lengths = length(length_vec);


temp = 300;
gnr_width = 7.5e-9;
rho_interlayer = 3e-1;
prob_backscattering = 0.0;
Ef = 0.2;
contact_resistance = 0;
defect_mfp = 1000e-9;

Rnet_mat = zeros(max_layers,num_lengths);
Reff_mat = zeros(max_layers,num_lengths);

for nind = 1:length(layer_vec)
    for lind = 1:num_lengths
        
        num_layers = layer_vec(nind);
        gnr_length = length_vec(lind);

        [Rnet Reff] = calc_gnr_params(num_layers,gnr_width,gnr_length,temp,defect_mfp,rho_interlayer,prob_backscattering,Ef,contact_resistance);
        Rnet_mat(nind,lind) = Rnet;
        Reff_mat(nind,lind) = Reff;
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
for pind = 1:length(layer_vec)
    plot(gnr_kumar_length*1e6,gnr_kumar_R(pind,:)/1e3)
end

set(gca,'ColorOrder',[0 0 0; 0 0 1; 0 1 0; 1 0 0]);
set(gca,'LineStyleOrder','--')
% plot(length_vec*1e6,Reff_mat/1e3)
for pind = 1:length(layer_vec)
    plot(length_vec*1e6,Reff_mat(pind,:)/1e3)
end
fixfigs(3,3,14,12)





