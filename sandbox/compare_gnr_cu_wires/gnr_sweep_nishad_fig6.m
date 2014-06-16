% Sweep mGNR params
clear all
close all

layer_vec = 1:1:10;
length_vec = (1:1:100)*1e-6;
num_lengths = length(length_vec);


temp = 300;
gnr_width = 7.8e-9;
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
R_com_top_mat = zeros(total_layers,num_lengths);
R_com_side_mat = zeros(total_layers,num_lengths);

delay_com_top_mat = zeros(total_layers,num_lengths);
delay_com_side_mat = zeros(total_layers,num_lengths);

C_com_mat = zeros(total_layers,num_lengths);
L_com_mat = zeros(total_layers,num_lengths);

height_dielectric = 16e-9;
epsrd = 2;

for nind = 1:length(layer_vec)
    for lind = 1:num_lengths
        
        num_layers = layer_vec(nind);
        gnr_length = length_vec(lind);

        [Rnet Reff] = calc_gnr_params(num_layers,gnr_width,gnr_length,temp,defect_mfp,rho_interlayer,prob_backscattering,Ef,contact_resistance);
        Rnet_mat(nind,lind) = Rnet;
        Reff_mat(nind,lind) = Reff;
        
        [delay_top delay_side R_top R_side L C] = calc_gnr_params_nishad(num_layers,gnr_width,gnr_length,temp,defect_mfp,rho_interlayer,prob_backscattering,Ef,contact_resistance,epsrd,height_dielectric);
        [delay_top_com delay_side_com R_top_com R_side_com L_com C_com] = calc_gnr_params_combined(num_layers,gnr_width,gnr_length,temp,defect_mfp,rho_interlayer,prob_backscattering,Ef,contact_resistance,epsrd,height_dielectric);
        
        R_nishad_top_mat(nind,lind) = R_top;
        R_nishad_side_mat(nind,lind) = R_side;
        
        R_com_top_mat(nind,lind) = R_top_com;
        R_com_side_mat(nind,lind) = R_side_com;
        
        C_com_mat(nind,lind) = C_com;
        L_com_mat(nind,lind) = L_com;
        
        delay_com_top_mat(nind,lind) = delay_top_com;
        delay_com_side_mat(nind,lind) = delay_side_com;
    end
end



%% Nishad fig 6: Delay


nind = 5;

figure(1)
clf
plot(length_vec*1e6,delay_com_top_mat(nind,:)*1e12,'b')
hold on
plot(length_vec*1e6,delay_com_side_mat(nind,:)*1e12,'r')
xlabel('Interconnect length ({\mu}m)')
ylabel('Delay (ps)')
set(gca,'yscale','log')
fixfigs(1,3,14,12)



%% Debug fig: Capacitance

figure(2)
clf
plot(length_vec*1e6, C_com_mat(nind,:)*1e12);
xlabel('Interconnect length ({\mu}m)')
ylabel('Capacitance (pF)')
set(gca,'yscale','log')
fixfigs(2,3,14,12)

figure(3)
clf
plot(length_vec*1e6, R_com_top_mat(nind,:)/1e3,'b');
hold on
plot(length_vec*1e6, R_com_side_mat(nind,:)/1e3,'r');
xlabel('Interconnect length ({\mu}m)')
ylabel('Resistance (k\Omega)')
%set(gca,'yscale','log')
fixfigs(3,3,14,12)

figure(4)
clf
plot(length_vec*1e6, L_com_mat(nind,:)*1e9,'b');
xlabel('Interconnect length ({\mu}m)')
ylabel('Inductance (nH)')
%set(gca,'yscale','log')
fixfigs(4,3,14,12)

%% Data from Kumar 2012 fig 14

k12_length_gp = 10:10:100;
delay_k12_tc0 = [ 2.99888	5.94557	9.58444	13.3884	17.5485	21.9287	26.5438	31.1234	35.917	41.4488 ];
delay_k12_sc0 = [ 2.99888	5.94557	9.4331	12.7641	16.4661	20.2512	24.1262	28.2887	32.6456	37.0787 ];
delay_k12_tc2 = [ 5.75929	13.177	21.9287	31.6228	42.7894	53.4701	66.8167	79.6016	93.3353	107.711 ];
delay_k12_sc2 = [ 4.39397	9.73821	16.2061	23.7452	32.1301	42.1138	53.4701	65.7616	78.3446	93.3353 ];
delay_k12_cu =  [ 5.31871	10.7141	16.4661	22.2806	29.2037	35.917	43.4759	50.9769	58.8282	67.8887 ];

figure(5)
clf
plot(k12_length_gp,delay_k12_tc0,'b')
hold on
plot(k12_length_gp,delay_k12_sc0,'r')
plot(k12_length_gp,delay_k12_tc2,'b--')
plot(k12_length_gp,delay_k12_sc2,'r--')
plot(k12_length_gp,delay_k12_cu, 'k')
set(gca,'yscale','log')
xlabel('length (GP)')
ylabel('delay (ps)')
fixfigs(5,3,14,12)






