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


%% Power vs scaling for different 3d configurations

figure(1)
clf
hold on
figure(2)
clf
hold on
colors = [ 0 0 0 ; 0 0 1; 0 1 0 ; 1 0 0 ];
linestyles = {'-', '--', ':'};
for nind = 1:num_stacks
    for bar_thick_ind = 1:num_barrier_thicknesses
        wire_flag_ind = 1;
        pow_tot = zeros(1,num_scaling_factors);
        pow_wire = zeros(1,num_scaling_factors);
        pow_rep = zeros(1,num_scaling_factors);

        pow_tot(1,:) = power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,:,bar_thick_ind,bar_res_ind);
        pow_wire(1,:) = wire_power(cind,dind,thind,nind,pind,:,wire_res_ind,wire_flag_ind,:,bar_thick_ind,bar_res_ind);
        pow_rep(1,:) = rep_power(cind,dind,thind,nind,pind,:,wire_res_ind,wire_flag_ind,:,bar_thick_ind,bar_res_ind);
        pow_comm = pow_wire + pow_rep;


        pow_logic = pow_tot - (pow_wire + pow_rep);
        pow_eff = pow_logic./pow_tot;
        pow_comm_ratio = pow_comm./pow_tot;
        pow_log_ratio = 1 - pow_comm_ratio;

        figure(1)
        plot(pow_comm_ratio,'color',colors(nind,:), 'linestyle', linestyles{bar_thick_ind})

        figure(2)
        plot(pow_tot, 'color', colors(nind,:), 'linestyle', linestyles{bar_thick_ind})
    end
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

%% Number of metal levels required for routing

figure(3)
clf
hold on
ymax = 0;
for bar_thick_ind = 1:num_barrier_thicknesses
    for nind = 1:num_stacks
        num_metal_levels = zeros(1, num_scaling_factors);

        for scaling_ind = 1:num_scaling_factors
            num_metal_levels(1,scaling_ind) = length(wire_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind}.pn);
        end

        plot(num_metal_levels, 'color', colors(nind,:), 'linestyle', linestyles{bar_thick_ind})
        ymax = max([ymax num_metal_levels]);
    end
end
%ylim([0 ymax+1])
set(gca,'Xtick',1:num_scaling_factors)
set(gca,'XtickLabel', {'22nm', '14nm', '10nm', '7nm', '5nm'} )
xlabel('Process Node')
ylabel('Number of Metal Levels')
fixfigs(3,2,14,12)


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
% figure(4)
% clf
% hold on
% for nind = 1:num_stacks
%     colors = [ 0 0 0 ; 0 0 1; 0 1 0 ; 1 0 0 ];
%     pow_tot = zeros(1,num_freqs);
%     pow_wire = zeros(1,num_freqs);
%     pow_rep = zeros(1,num_freqs);
%     
%     pow_tot(1,:) = power(cind,dind,thind,nind,pind,:,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind);
%     pow_wire(1,:) = wire_power(cind,dind,thind,nind,pind,:,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind);
%     pow_rep(1,:) = rep_power(cind,dind,thind,nind,pind,:,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind);
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
%     fixfigs(4,2,14,12)
% end
%     
%     
%     
% power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind) = core.power.total;
% power_density(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind) = core.power.density;
% 
% wire_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind) = core.power.wiring;
% rep_power(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind) = core.power.repeater;


%%  Max frequency vs ILD
% Sweeping number of tiers, epsrd, and using thermall-limited frequency
% search
% 
% figure(5)
% clf
% hold on
% 
% for nind = 1:num_stacks
%     colors = [ 0 0 0 ; 0 0 1; 0 1 0 ; 1 0 0 ];
%     fr_vec1 = zeros(1,num_perms);
%     fr_vec2 = zeros(1,num_perms);
%     cind = 1;
%     fr_vec1(1,:) = freq(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind)/1e9;
%     cind = 2;
%     fr_vec2(1,:) = freq(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind)/1e9;
%     
%     plot(rel_permittivities,fr_vec1,'color',colors(nind,:),'linestyle','-')
%     plot(rel_permittivities,fr_vec2,'color',colors(nind,:),'linestyle','--')
% end
% xlabel('ILD Relative Permittivity')
% ylabel('Maximum Frequency')
% fixfigs(5,2,14,12)
% 
% figure(6)
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
%     pow_vec(1,:) = power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind);
%     comm_pow_vec(1:end) = wire_power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind) + rep_power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind);
%     p_ave_vec(nind) = mean(pow_vec);
%     comm_pow_ave_vec(nind) = mean(comm_pow_vec);
%     
%     cind = 2;
%     pow_vec2 = zeros(1,num_perms);
%     comm_pow_vec2 = zeros(1,num_perms);
%     pow_vec2(1,:) = power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind);
%     comm_pow_vec2(1:end) = wire_power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind) + rep_power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind);
%     p_ave_vec2(nind) = mean(pow_vec2);
%     comm_pow_ave_vec2(nind) = mean(comm_pow_vec2);
%     
%     plot(rel_permittivities,pow_vec,'color',colors(nind,:),'linestyle','-')
%     plot(rel_permittivities,pow_vec2,'color',colors(nind,:),'linestyle','--')
% end
% xlabel('ILD Relative Permittivity')
% ylabel('Power (W)')
% fixfigs(6,2,14,12)
% 
% 
% pfrac_vec = comm_pow_ave_vec./p_ave_vec;
% pfrac_vec2 = comm_pow_ave_vec2./p_ave_vec2;
% 
% figure(7)
% clf
% hold on
% plot(tiers,p_ave_vec,'r')
% plot(tiers,p_ave_vec2,'b')
% xlabel('Tiers')
% ylabel('Power (W)')
% fixfigs(7,2,14,12)
% 
% 
% pmat = [p_ave_vec ; p_ave_vec2]';
% 
% figure(8)
% clf
% b = bar(pmat,1,'grouped');
% colormap jet
% set(gca,'xticklabel',{'1','2','4','8'})
% xlabel('Number of tiers')
% ylabel('Power (W)')
% b(1).FaceColor = 'blue';
% b(2).FaceColor = 'yellow';
% fixfigs(8,2,14,12)
% 
% 
% figure(9)
% clf
% pmat = [pfrac_vec ; pfrac_vec2]';
% b = bar(pmat,1,'grouped');
% colormap jet
% set(gca,'xticklabel',{'1','2','4','8'})
% xlabel('Number of tiers')
% ylabel('Comm Power Fraction')
% b(1).FaceColor = 'blue';
% b(2).FaceColor = 'yellow';
% fixfigs(9,2,14,12)
% 
% 
% % Plots
% figure(10)
% clf
% hold on
% for nind = 1:num_stacks
%     colors = [ 0 0 0 ; 0 0 1; 0 1 0 ; 1 0 0 ];
%     fr_vec1 = zeros(1,num_perms);
%     fr_vec2 = zeros(1,num_perms);
%     pow_vec = zeros(1,num_perms);
%     pow_vec2 = zeros(1,num_perms);
%     cind = 1;
%     fr_vec1(1,:) = freq(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind);
%     pow_vec(1,:) = power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind);
%     epc_vec = pow_vec./fr_vec1;
%     
%     cind = 2;
%     fr_vec2(1,:) = freq(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind);
%     pow_vec2(1,:) = power(cind,dind,thind,nind,:,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind);
%     epc_vec2 = pow_vec2./fr_vec2;
%     
%     plot(rel_permittivities,epc_vec*1e9,'color',colors(nind,:),'linestyle','-')
%     plot(rel_permittivities,epc_vec2*1e9,'color',colors(nind,:),'linestyle','--')
% end
% xlabel('ILD Relative Permittivity')
% ylabel('Energy per cycle (nJ)')
% fixfigs(10,2,14,12)

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
% pvec(1,:) = power(cind,dind,thind,:,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind);
% pdens_vec(1,:) = power_density(cind,dind,thind,:,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind)/100^2;
% 
% figure(11)
% clf
% [ax, h1, h2] = plotyy(tiers,pvec,tiers,pdens_vec);
% set(ax(1),'ycolor','k')
% set(ax(2),'ycolor','k')
% ax(1).FontSize = 12;
% ax(2).FontSize = 12;
% %ax(1).YLabel.String = 'Power (W)';
% %ax(2).YLabel.String = 'Power Density (W/cm^2)';
% xlabel('Number of tiers')
% fixfigs(11,2,14,12)

%% Max frequency vs barrier thickness across nodes

figure(12)
clf
hold on

colors = [ 0 0 0 ; 0 0 1; 0 1 0 ; 1 0 0 ];
linestyles = {'-','--',':'};

for nind = 1:num_stacks
    for bar_thick_ind = 1:num_barrier_thicknesses
        fr_vec = zeros(1,num_scaling_factors);
        fr_vec(1,:) = freq(cind,dind,thind,nind,cind,freq_ind,wire_res_ind,wire_flag_ind,:,bar_thick_ind,bar_res_ind)/1e9;
        plot(fr_vec,'color',colors(nind,:),'linestyle',linestyles{bar_thick_ind})
    end
end
set(gca,'Xtick',1:num_scaling_factors)
set(gca,'XtickLabel', {'22nm', '14nm', '10nm', '7nm', '5nm'} )
xlabel('Process Node')
ylabel('Maximum Frequency')
fixfigs(12,2,14,12)


%% Leakage power for different scaling scenarios

figure(13)
clf
hold on
figure(14)
clf
hold on
colors = [ 0 0 0 ; 0 0 1; 0 1 0 ; 1 0 0 ];
linestyles = {'-', '--', ':'};
for nind = 1:num_stacks
    for bar_thick_ind = 1:num_barrier_thicknesses
        pow_leakage = zeros(1,num_scaling_factors);
        pow_leakage(1,:) = leakage_power(cind,dind,thind,nind,pind,:,wire_res_ind,wire_flag_ind,:,bar_thick_ind,bar_res_ind);
        
        pow_dynamic = zeros(1,num_scaling_factors);
        pow_dynamic(1,:) = dynamic_power(cind,dind,thind,nind,pind,:,wire_res_ind,wire_flag_ind,:,bar_thick_ind,bar_res_ind);

        figure(13)
        plot(pow_dynamic,'color',colors(nind,:), 'linestyle', linestyles{bar_thick_ind})

        figure(14)
        plot(pow_leakage, 'color', colors(nind,:), 'linestyle', linestyles{bar_thick_ind})
    end
end
figure(13)
set(gca,'Xtick',1:num_scaling_factors)
set(gca,'XtickLabel', {'22nm', '14nm', '10nm', '7nm', '5nm'} )
xlabel('Process Node')
ylabel('Dynamic Power (W)')
fixfigs(13,2,14,12)

figure(14)
set(gca,'Xtick',1:num_scaling_factors)
set(gca,'XtickLabel', {'22nm', '14nm', '10nm', '7nm', '5nm'} )
xlabel('Process Node')
ylabel('Leakage Power (W)')
fixfigs(14,2,14,12)


%% Number of metal levels vs wire resistivity and scaling

% for nind = 1:num_stacks
%     fignum = 15 + nind -1;
%     figure(fignum)
%     clf
%     hold on
%     
%     num_metal_levels = zeros(num_wire_resistivities, num_scaling_factors);
%     for scaling_ind = 1:num_scaling_factors
%         for wire_res_ind = 1:num_wire_resistivities
%             num_metal_levels(wire_res_ind,scaling_ind) = length(wire_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind}.pn);
%         end
%     end
% 
%     surf( 1:num_scaling_factors, wire_resistivities*1e9, num_metal_levels);
%     view(2)
%     set(gca,'Xtick',1:num_scaling_factors)
%     set(gca,'XtickLabel', {'22nm', '14nm', '10nm', '7nm', '5nm'} )
%     xlabel('Process Node')
%     ylabel('Wire Resistivity (Ohm-nm)')
% 
%     colormap jet
%     colorbar
%     fixfigs(fignum,3,14,12)
% end


