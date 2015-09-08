function plot_sweep_data( sweep, sweep_data, simulation)

%% Unpack sweep data
power = sweep_data.power;
power_density = sweep_data.power_density;
freq = sweep_data.freq;

wire_power = sweep_data.wire_power;
rep_power = sweep_data.rep_power;
temp = sweep_data.temp;
thickness = sweep_data.thickness;
npads = sweep_data.npads;


ild_cell = sweep_data.ild_cell;
psn_cell = sweep_data.psn_cell;
power_cell = sweep_data.power_cell;
wire_cell = sweep_data.wire_cell;
chip_cell = sweep_data.chip_cell;
tsv_cell = sweep_data.tsv_cell;

Ltsv_m2 = sweep_data.Ltsv_m2;
cap_density = sweep_data.cap_density;

num_stacks = length(sweep.tiers);
num_perms = length(sweep.rel_permittivities);
num_thicks = length(sweep.thicknesses);

if (simulation.freq_binsearch == 1)
    num_freqs = 1; % skip freq loops
else
    num_freqs = length(sweep.frequencies);
end
num_cooling_configs = length(sweep.heat_fluxes);
num_decaps = length(sweep.decap_ratios);
num_wire_resistivities = length(sweep.wire_resistivities);
num_wire_flags = length(sweep.wire_material_flags);
num_scaling_factors = length(sweep.scaling_factors);

cind = num_cooling_configs;
dind = num_decaps;
thind = num_thicks;
nind = num_stacks;
pind = num_perms;
freq_ind = num_freqs;
wire_res_ind = num_wire_resistivities;
wire_flag_ind = num_wire_flags;
scaling_ind = num_scaling_factors;

tiers = sweep.tiers;
thicknesses = sweep.thicknesses;
force_thickness = sweep.force_thickness;
rel_permittivities = sweep.rel_permittivities;
frequencies = sweep.frequencies;
heat_fluxes = sweep.heat_fluxes;
decap_ratios = sweep.decap_ratios;
wire_resistivities = sweep.wire_resistivities;
wire_material_flags = sweep.wire_material_flags;
scaling_factors = sweep.scaling_factors;


%% Power vs scaling for different 3d configurations

figure(1)
clf
hold on
figure(2)
clf
hold on
for nind = 1:num_stacks
    colors = [ 0 0 0 ; 0 0 1; 0 1 0 ; 1 0 0 ];
    
    wire_flag_ind = 1;
    pow_tot = zeros(1,num_scaling_factors);
    pow_wire = zeros(1,num_scaling_factors);
    pow_rep = zeros(1,num_scaling_factors);
    
    pow_tot(1,:) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,:);
    pow_wire(1,:) = wire_power(cind,dind,thind,nind,pind,:,wire_res_ind,wire_flag_ind,:);
    pow_rep(1,:) = rep_power(cind,dind,thind,nind,pind,:,wire_res_ind,wire_flag_ind,:);
    pow_comm = pow_wire + pow_rep;
    
    
    pow_logic = pow_tot - (pow_wire + pow_rep);
    pow_eff = pow_logic./pow_tot;
    pow_comm_ratio = pow_comm./pow_tot;
    pow_log_ratio = 1 - pow_comm_ratio;
    
    figure(1)
    plot(pow_comm_ratio,'color',colors(nind,:))
    
    figure(2)
    plot(pow_tot, 'color', colors(nind,:) )
    
%     wire_flag_ind = 2;
%     pow_tot = zeros(1,num_scaling_factors);
%     pow_wire = zeros(1,num_scaling_factors);
%     pow_rep = zeros(1,num_scaling_factors);
%     
%     pow_tot(1,:) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,:);
%     pow_wire(1,:) = wire_power(cind,dind,thind,nind,pind,:,wire_res_ind,wire_flag_ind,:);
%     pow_rep(1,:) = rep_power(cind,dind,thind,nind,pind,:,wire_res_ind,wire_flag_ind,:);
%     pow_comm = pow_wire + pow_rep;
%     
%     
%     pow_logic = pow_tot - (pow_wire + pow_rep);
%     pow_eff = pow_logic./pow_tot;
%     pow_comm_ratio = pow_comm./pow_tot;
%     pow_log_ratio = 1 - pow_comm_ratio;
%     
%     figure(1)
%     plot(pow_comm_ratio,'color',colors(nind,:),'linestyle','--')
%     
%     figure(2)
%     plot(pow_tot, 'color', colors(nind,:),'linestyle','--' )
end
figure(1)
set(gca,'Xtick',1:num_scaling_factors)
set(gca,'XtickLabel', {'22nm', '14nm', '10nm', '7nm', '5nm'} )
xlabel('Process Node')
ylabel('On-chip Communication Power Fraction')
fixfigs(1,2,14,12)

figure(2)
set(gca,'Xtick',1:num_scaling_factors)
set(gca,'XtickLabel', {'22nm', '14nm', '10nm', '7nm', '5nm'} )
xlabel('Process Node')
ylabel('Power Consumption (W)')
fixfigs(2,2,14,12)


%% Power efficiency vs frequency
% Inputs
% tiers = [1 2 4 8];
% thicknesses = [10e-6];
% force_thickness = 1;
% rel_permittivities = 3.0;
% frequencies = linspace(0.1e9, 5e9, 1e2);
% heat_fluxes = [ h_air ];
% decap_ratios = [0.1];%[0.01 0.1 1];
% %wire_resistivities = [rho_ag rho_cu rho_au rho_al rho_w rho_ni];
% wire_resistivities = [rho_cu];
% wire_material_flags = {'00'}; % binary strings. bit1 = use_graphene, bit0 = use alt_em_mat
% scaling_factor = [1];
% 
% Plots
% figure(1)
% clf
% hold on
% for nind = 1:num_stacks
%     colors = [ 0 0 0 ; 0 0 1; 0 1 0 ; 1 0 0 ];
%     pow_tot = zeros(1,num_freqs);
%     pow_wire = zeros(1,num_freqs);
%     pow_rep = zeros(1,num_freqs);
%     
%     pow_tot(1,:) = power(cind,dind,thind,nind,pind,:,wire_res_ind,wire_flag_ind,scaling_ind);
%     pow_wire(1,:) = wire_power(cind,dind,thind,nind,pind,:,wire_res_ind,wire_flag_ind,scaling_ind);
%     pow_rep(1,:) = rep_power(cind,dind,thind,nind,pind,:,wire_res_ind,wire_flag_ind,scaling_ind);
%     pow_comm = pow_wire + pow_rep;
%     
%     
%     pow_logic = pow_tot - (pow_wire + pow_rep);
%     pow_eff = pow_logic./pow_tot;
%     pow_comm_ratio = pow_comm./pow_tot;
%     pow_log_ratio = 1 - pow_comm_ratio;
%     
%     plot(frequencies/1e9,pow_comm_ratio,'color',colors(nind,:))
%     xlim([0.5 5])
%     xlabel('Clock Frequency (GHz)')
%     ylabel('On-chip Communication Power Fraction')
%     %set(gca,'xscale','log')
%     fixfigs(1,2,14,12)
% end
%     
%     
%     
% power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.power.total;
% power_density(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.power.density;
% 
% wire_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.power.wiring;
% rep_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) = core.power.repeater;


%%  Max frequency vs ILD
% % Inputs
% tiers = [1 2 4 8];
% thicknesses = [10e-6];
% force_thickness = 1;
% rel_permittivities = linspace(1,4,41);
% frequencies = [3e8 1e10]; % if simulation.freq_binsearch is set, (1) is min freq and (2) is max freq
% heat_fluxes = [ h_air h_water];
% decap_ratios = [0.1];%[0.01 0.1 1];
% %wire_resistivities = [rho_ag rho_cu rho_au rho_al rho_w rho_ni];
% wire_resistivities = [rho_cu];
% wire_material_flags = {'00'}; % binary strings. bit1 = use_graphene, bit0 = use alt_em_mat
% scaling_factor = [1];
% 
% figure(1)
% clf
% hold on
% 
% for nind = 1:num_stacks
%     colors = [ 0 0 0 ; 0 0 1; 0 1 0 ; 1 0 0 ];
%     fr_vec1 = zeros(1,num_perms);
%     fr_vec2 = zeros(1,num_perms);
%     cind = 1;
%     fr_vec1(1,:) = freq(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind)/1e9;
%     cind = 2;
%     fr_vec2(1,:) = freq(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind)/1e9;
%     
%     plot(rel_permittivities,fr_vec1,'color',colors(nind,:),'linestyle','-')
%     plot(rel_permittivities,fr_vec2,'color',colors(nind,:),'linestyle','--')
% end
% xlabel('ILD Relative Permittivity')
% ylabel('Maximum Frequency')
% fixfigs(1,2,14,12)
% 
% figure(2)
% clf
% hold on
% p_ave_vec = zeros(1,num_stacks);
% p_ave_vec2 = zeros(1,num_stacks);
% comm_pow_ave_vec = zeros(1,num_stacks);
% comm_pow_ave_vec2 = zeros(1,num_stacks);
% for nind = 1:num_stacks
%     colors = [ 0 0 0 ; 0 0 1; 0 1 0 ; 1 0 0 ];
%     cind = 1;
%     pow_vec = zeros(1,num_perms);
%     comm_pow_vec = zeros(1,num_perms);
%     pow_vec(1,:) = power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind);
%     comm_pow_vec(1:end) = wire_power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) + rep_power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind);
%     p_ave_vec(nind) = mean(pow_vec);
%     comm_pow_ave_vec(nind) = mean(comm_pow_vec);
%     
%     cind = 2;
%     pow_vec2 = zeros(1,num_perms);
%     comm_pow_vec2 = zeros(1,num_perms);
%     pow_vec2(1,:) = power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind);
%     comm_pow_vec2(1:end) = wire_power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind) + rep_power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind);
%     p_ave_vec2(nind) = mean(pow_vec2);
%     comm_pow_ave_vec2(nind) = mean(comm_pow_vec2);
%     
%     plot(rel_permittivities,pow_vec,'color',colors(nind,:),'linestyle','-')
%     plot(rel_permittivities,pow_vec2,'color',colors(nind,:),'linestyle','--')
% end
% xlabel('ILD Relative Permittivity')
% ylabel('Power (W)')
% fixfigs(2,2,14,12)
% 
% 
% pfrac_vec = comm_pow_ave_vec./p_ave_vec;
% pfrac_vec2 = comm_pow_ave_vec2./p_ave_vec2;
% 
% figure(3)
% clf
% hold on
% plot(tiers,p_ave_vec,'r')
% plot(tiers,p_ave_vec2,'b')
% xlabel('Tiers')
% ylabel('Power (W)')
% fixfigs(3,2,14,12)
% 
% 
% pmat = [p_ave_vec ; p_ave_vec2]';
% 
% figure(4)
% clf
% b = bar(pmat,1,'grouped');
% colormap jet
% set(gca,'xticklabel',{'1','2','4','8'})
% xlabel('Number of tiers')
% ylabel('Power (W)')
% b(1).FaceColor = 'blue';
% b(2).FaceColor = 'yellow';
% fixfigs(4,2,14,12)
% 
% 
% figure(5)
% clf
% pmat = [pfrac_vec ; pfrac_vec2]';
% b = bar(pmat,1,'grouped');
% colormap jet
% set(gca,'xticklabel',{'1','2','4','8'})
% xlabel('Number of tiers')
% ylabel('Comm Power Fraction')
% b(1).FaceColor = 'blue';
% b(2).FaceColor = 'yellow';
% fixfigs(5,2,14,12)
% 
% 
% % Plots
% figure(6)
% clf
% hold on
% for nind = 1:num_stacks
%     colors = [ 0 0 0 ; 0 0 1; 0 1 0 ; 1 0 0 ];
%     fr_vec1 = zeros(1,num_perms);
%     fr_vec2 = zeros(1,num_perms);
%     pow_vec = zeros(1,num_perms);
%     pow_vec2 = zeros(1,num_perms);
%     cind = 1;
%     fr_vec1(1,:) = freq(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind);
%     pow_vec(1,:) = power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind);
%     epc_vec = pow_vec./fr_vec1;
%     
%     cind = 2;
%     fr_vec2(1,:) = freq(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind);
%     pow_vec2(1,:) = power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind);
%     epc_vec2 = pow_vec2./fr_vec2;
%     
%     plot(rel_permittivities,epc_vec*1e9,'color',colors(nind,:),'linestyle','-')
%     plot(rel_permittivities,epc_vec2*1e9,'color',colors(nind,:),'linestyle','--')
% end
% xlabel('ILD Relative Permittivity')
% ylabel('Energy per cycle (nJ)')
% fixfigs(6,2,14,12)

% epc_mat = zeros(4,num_stacks);

%% Power consumption vs power density

% % Inputs
% tiers = 1:8;
% thicknesses = [10e-6];
% force_thickness = 1;
% rel_permittivities = [3.0];
% frequencies = fmax_core; % if simulation.freq_binsearch is set, (1) is min freq and (2) is max freq
% heat_fluxes = [ h_air];
% decap_ratios = [0.1];%[0.01 0.1 1];
% %wire_resistivities = [rho_ag rho_cu rho_au rho_al rho_w rho_ni];
% wire_resistivities = [rho_cu];
% wire_material_flags = {'00'}; % binary strings. bit1 = use_graphene, bit0 = use alt_em_mat
% scaling_factor = [1];
% 
% 
% pvec = zeros(1,num_stacks);
% pdens_vec = zeros(1,num_stacks);
% pvec(1,:) = power(cind,dind,thind,:,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind);
% pdens_vec(1,:) = power_density(cind,dind,thind,:,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind)/100^2;
% 
% figure(1)
% clf
% [ax, h1, h2] = plotyy(tiers,pvec,tiers,pdens_vec);
% set(ax(1),'ycolor','k')
% set(ax(2),'ycolor','k')
% ax(1).FontSize = 12;
% ax(2).FontSize = 12;
% %ax(1).YLabel.String = 'Power (W)';
% %ax(2).YLabel.String = 'Power Density (W/cm^2)';
% xlabel('Number of tiers')
% fixfigs(1,2,14,12)
