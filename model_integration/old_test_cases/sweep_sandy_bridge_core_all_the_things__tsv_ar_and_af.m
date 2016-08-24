%% Simulation parameters
simulation.use_joyner = 0;
simulation.redo_wiring_after_repeaters = 0;
simulation.topdown_WLARI = 1; % Use topdown simultaneous WLA and RI (0 = use standard bottom-up optimal WLA, followed by one pass of RI)
simulation.skip_psn_loops = 0; % Skip PSN TSV homing for faster debug
simulation.draw_thermal_map = 0; % Plot thermal profile of each chip
simulation.print_thermal_data = 0; % Output max temp in each layer to console
simulation.separate_wiring_tiers = 1; % 1 = Each logic plane will have its own wiring tiers between it and the next logic plane
                                      % 0 = All metal layers for entire device will be routed on top of entire 3D stack

%% Logic core parameters
Ng_core = 86e6/4;
Ach_mm2_core = 18.5;
gate_pitch_core = 465e-9*2;
min_pitch_core = 112.5e-9;
fmax_core = 3.6e9;
w_trans = 32e-9;
Vdd_core = 1.25;

%% Thermal parameters
% %the heat transfer coefficient
% % r = 1/(hA); A is the size of top surface area
% % the cooling capability of the top heatsink; 20000, 1cm*1cm, means:
% % 0.5 W/K
% % h = q/dT - q = heat flux (W/m^2)
% heat.up = 20000;
% 
% % Bottom surface heat transfer coefficient
% % This parameter controls the area directly BELOW the bottom chip
% % If the interposer is larger than the bottom chip, heat.d controls the
% % rest of the area
% % Microfluidic heat sinks are assumed to be as large as the chip in the interposer
% heat.down = 5;  
% 
% % Heat transfer coefficient for the interposer area NOT directly underneath
% % the chip(s)
% heat.d = 5;
% 
% % Side surface heat coefficient, usually near adiabatic
% heat.side = 5;
% 
% heat.Ta = 298; % ambient temperature


r_air = 1/1.825; %K/W for a 1cm^2 HS
r_water = 1/4.63; %K/W for a 1cm^2 HS
A_hs = (1e-2)^2; % 1 cm^2

h_air = 1/(r_air*A_hs);
h_water = 1/(r_water*A_hs);
h_package = 5; % it sucks

heat.up = h_air;        % above chip
heat.down = h_package;     % directly beneath chip
heat.d = h_package;        % package, not under chip
heat.side = h_package;          % side
heat.Ta = 298; % ambient temperature

heat.interposer_thickness = 200e-6; % (m) Thickness of the interposer below the 3D stack
heat.bump_thickness = 40e-6;        % (m) Microbump thickness (between interposer and bottom chip of 3D stack)
heat.underfill_thickness = 1e-6;    % (m) Thickness of underfill material between each die in the 3D stack
heat.tim_thickness = 5e-6;          % (m) Thickness of thermal interface material between top chip in stack and heat sink

heat.material_IDs = [ 2 9 3];



%%
rent_exp_logic = 0.6;
rent_exp_mem = 0.4;
rent_exp_gpu = 0.55;

%% 
tiers = 1:8;
thicknesses = [50e-6];
force_thickness = 0;
rel_permittivities = [3.9];
frequencies = [fmax_core];
heat_fluxes = [ h_air ];
%decap_ratios = [1e-2 1e-1:1e-1:1];
%decap_ratios = [0.01 0.1 1.0];
decap_ratios = 0.1;


num_stacks = length(tiers);
num_perms = length(rel_permittivities);
num_thicks = length(thicknesses);
num_freqs = length(frequencies);
num_cooling_configs = length(heat_fluxes);
num_decaps = length(decap_ratios);
total_configs = num_stacks * num_perms * num_thicks * num_freqs * num_cooling_configs * num_decaps;

power = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs);
power_density = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs);

wire_power = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs);
rep_power = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs);
temp = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs);
thickness = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs);
npads = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs);
cap_density = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs);
Ltsv_m2 = zeros(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs);


ild_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs);
psn_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs);
power_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs);
wire_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs);
chip_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs);
tsv_cell = cell(num_cooling_configs,num_decaps,num_thicks,num_stacks,num_perms,num_freqs);

t_sweep_start = cputime;
cur_config = 0;
for cind = 1:num_cooling_configs
    for dind = 1:num_decaps
        for thind = 1:num_thicks
            for nind = 1:length(tiers)
                for pind = 1:num_perms
                    for freq_ind = 1:num_freqs
                        cur_config = cur_config + 1;
                        die_thickness = thicknesses(thind);
                        num_layers_per_block = tiers(nind);
                        epsrd = rel_permittivities(pind);
                        fmax_core = frequencies(freq_ind);

                        fprintf('===============================\n')
                        fprintf('==   cooling: %d/%d \t=====\n',cind,num_cooling_configs);
                        fprintf('==     decap: %d/%d \t=====\n',dind,num_decaps);
                        fprintf('== thickness: %d/%d \t=====\n',thind,num_thicks);
                        fprintf('==     Tiers: %d/%d \t=====\n',nind,num_stacks);
                        fprintf('==     epsrd: %d/%d \t=====\n',pind,num_perms);
                        fprintf('==      freq: %d/%d \t=====\n',freq_ind,num_freqs);
                        fprintf('==   Overall: %d/%d \t=====\n',cur_config,total_configs);
                        fprintf('===============================\n')

                        %% define parameters

                        [core.chip core.transistor core.gate core.tsv core.wire core.psn] = generate_basic_processor_settings(rent_exp_logic,num_layers_per_block,Ng_core,Ach_mm2_core,gate_pitch_core,min_pitch_core,Vdd_core,fmax_core,w_trans);
                        %core.psn.mismatch_tolerance = 0.01;
                        %% Tweak wiring parameters
                        core.wire.repeater_fraction = [0.3]; % 1 is default from gen_basic_proc_settings
                        core.wire.routing_efficiency = [0.6]; % 0.4 is default from gen_basic_proc_settings
                        core.wire.use_graphene = 0;
                        simulation.force_thickness = force_thickness;
                        core.chip.thickness_nominal = die_thickness;
                        core.wire.dielectric_epsr = epsrd;
                        core.psn.decap_area_fraction = decap_ratios(dind);
                        
                        if (die_thickness < 30e-6) % for monolithic-scale chips use thin SiO2 layer rather than underfill
                            heat.interposer_thickness = 200e-6; % (m) Thickness of the interposer below the 3D stack
                            heat.bump_thickness = 40e-6;        % (m) Microbump thickness (between interposer and bottom chip of 3D stack)
                            heat.underfill_thickness = 0.2e-6;    % (m) Thickness of underfill material between each die in the 3D stack
                            heat.tim_thickness = 5e-6;          % (m) Thickness of thermal interface material between top chip in stack and heat sink
                            heat.material_IDs = [ 2 9 5];
                        else % for standard die stacking go ahead and use regular underfill
                            heat.interposer_thickness = 200e-6; % (m) Thickness of the interposer below the 3D stack
                            heat.bump_thickness = 40e-6;        % (m) Microbump thickness (between interposer and bottom chip of 3D stack)
                            heat.underfill_thickness = 5e-6;    % (m) Thickness of underfill material between each die in the 3D stack
                            heat.tim_thickness = 5e-6;          % (m) Thickness of thermal interface material between top chip in stack and heat sink
                            heat.material_IDs = [ 2 9 3];
                        end

                        heat.up = heat_fluxes(cind);        % above chip

                        %% calculate block parameters
                        [core.chip core.power core.tsv core.wire core.repeater core.psn] = codesign_block(core.chip,core.tsv,core.gate,core.transistor,core.wire,heat,core.psn,simulation);

                        power(cind,dind,thind,nind,pind,freq_ind) = core.power.total;
                        power_density(cind,dind,thind,nind,pind,freq_ind) = core.power.density;

                        wire_power(cind,dind,thind,nind,pind,freq_ind) = core.power.wiring;
                        rep_power(cind,dind,thind,nind,pind,freq_ind) = core.power.repeater;
                        temp(cind,dind,thind,nind,pind,freq_ind) = core.chip.temperature;
                        thickness(cind,dind,thind,nind,pind,freq_ind) = core.chip.thickness;
                        npads(cind,dind,thind,nind,pind,freq_ind) = core.psn.Npads;


                        ild_cell{cind,dind,thind,nind,pind,freq_ind} = core.chip.iidf;
                        psn_cell{cind,dind,thind,nind,pind,freq_ind} = core.psn;
                        power_cell{cind,dind,thind,nind,pind,freq_ind} = core.power;
                        wire_cell{cind,dind,thind,nind,pind,freq_ind} = core.wire;
                        chip_cell{cind,dind,thind,nind,pind,freq_ind} = core.chip;
                        tsv_cell{cind,dind,thind,nind,pind,freq_ind} = core.tsv;

                        Ltsv_m2(cind,dind,thind,nind,pind,freq_ind) = psn_cell{cind,dind,thind,nind,pind,freq_ind}.Ltsv/psn_cell{cind,dind,thind,nind,pind,freq_ind}.l_unit_cell^2;
                        cap_density(cind,dind,thind,nind,pind,freq_ind) = psn_cell{cind,dind,thind,nind,pind,freq_ind}.cap_density;
                    end
                end
            end
        end
    end
end

t_sweep_stop = cputime;
fprintf('\nTotal time elapsed for parameter sweep: %.3g seconds\n\n',(t_sweep_stop-t_sweep_start));




%% Power consumption vs tier number
%npads(cind,dind,thind,nind,pind,freq_ind)

% f1 = figure(1);
% clf
% 
% Yplot = [];
% for nind = 1:num_stacks
%     pvec = zeros(1,num_perms);
%     pvec(1:end) = power(1,1, 1 ,nind, : ,1);
%     pvec = fliplr(pvec);
%     Yplot = [Yplot; pvec];
%     
% end
% b = bar(Yplot,1.0,'group');
% set(b(1),'FaceColor',[0.00 0.00 0.85])
% set(b(2),'FaceColor',[0.00 0.85 0.00])
% set(b(3),'FaceColor',[0.90 0.90 0.00])
% set(b(4),'FaceColor',[0.85 0.00 0.00])
% % pow1 = zeros(1,length(tiers));
% % pow2 = pow1;
% % pow1(1:end) = power(1,1, 1 ,:, end ,1);
% % pow2(1:end) = power(1,1, 1 ,:, 1 ,1);
% % plot(tiers,pow1,'b');
% % hold on
% % plot(tiers,pow2,'r');
% xlabel('Number of tiers')
% ylabel('Total power (W)')
% xlim([0 9])
% fixfigs(1,3,14,12)
% legend({'\epsilon_r = 3.9','\epsilon_r = 3.0','\epsilon_r = 2.0','\epsilon_r = 1.0'},'fontweight','bold')
% legend('boxoff')
% 
% 
% f2 = figure(2);
% clf
% hold on
% col = [0 0 0; 0 0 1; 0 1 0 ; 1 0 0];
% for pind = num_perms:-1:1
%     pvec = zeros(1,num_stacks);
%     pvec(1:end) = power(1,1, 1 ,:, pind ,1);
%     %pvec = fliplr(pvec);
%     plot(tiers,pvec,'linestyle','-','color',col((num_perms-pind+1),:))
%     
% end
% ylim([0 35])
% xlabel('Number of tiers')
% ylabel('Total power (W)')
% fixfigs(2,3,14,12)

% tsv_ars = [ 5 10 20];
% tsv_area_fracs = [ 0.01 0.1 0.1];
% heat_fluxes = tsv_ars;
% frequencies = tsv_area_fracs;

%npads(cind,dind,thind,nind,pind,freq_ind)
%ar = h/c
%af = f
f3 = figure(3);
clf
hold on
col = [0 0 0; 0 0 1; 0 1 0 ; 1 0 0];
for cind = 1:num_cooling_configs
    pvec = zeros(1,num_stacks);
    pvec(1:end) = thickness(cind,1, 1 ,:, 1 ,1)*1e6;
    %pvec = fliplr(pvec);
    plot(tiers(2:end),pvec(2:end),'linestyle','-','color',col(cind,:))
    
end
%ylim([0 35])
xlabel('Number of tiers')
ylabel('Maximum die thickness (microns)')
fixfigs(3,3,14,12)


%npads(cind,dind,thind,nind,pind,freq_ind)
%ar = h/c
%af = f
f4 = figure(4);
clf
hold on
col = [0 0 0; 0 0 1; 0 1 0 ; 1 0 0];
sty = {'-' '-' '--'};
for freq_ind = 1:num_freqs
    pvec = zeros(1,num_stacks);
    pvec(1:end) = thickness(2,1, 1 ,:, 1 ,freq_ind)*1e6;
    %pvec = fliplr(pvec);
    plot(tiers(2:end),pvec(2:end),'linestyle','-','color',col(freq_ind,:))
    
end
%ylim([0 35])
xlabel('Number of tiers')
ylabel('Maximum die thickness (microns)')
fixfigs(4,3,14,12)

f5 = figure(5);
clf
Yplot = [];
for nind = [2 4 6 8]
    pvec = zeros(1,num_cooling_configs*2);
    pvec(1:num_cooling_configs) = thickness(:,1, 1 ,nind, 1 ,1)*1e6;
    pvec(num_cooling_configs+1:end) = thickness(:,1, 1 ,nind, 1 ,2)*1e6;
    pvec = fliplr(pvec);
    Yplot = [Yplot; pvec];
    
end
b = bar(Yplot,1.0,'group');
xlabel('Number of Tiers')
ylabel('Maximum Die Thickness (microns)')
set(gca,'yscale','log')
fixfigs(5,3,14,12)
%legend({'\epsilon_r = 3.9','\epsilon_r = 3.0','\epsilon_r = 2.0','\epsilon_r = 1.0'},'fontweight','bold')
%legend('boxoff')
set(gca,'xtickmode','manual')
set(gca,'xtick',1:4)
set(gca,'xticklabel',2:2:8)




%%
% %% Plot number of PSN TSVs for different configurations
% 
% %npads(cind,dind,thind,nind,pind,freq_ind)
% Yplot = [ npads(1,2, 1 ,1, 1 ,1)   npads(1,2, 2 ,1, 1 ,1)   npads(1,2, 3 ,1, 1 ,1)   npads(1,2, 1 ,1, 4 ,1)   npads(1,2, 2 ,1, 4 ,1)   npads(1,2, 3 ,1, 4 ,1);
%           npads(1,2, 1 ,2, 1 ,1)   npads(1,2, 2 ,2, 1 ,1)   npads(1,2, 3 ,2, 1 ,1)   npads(1,2, 1 ,2, 4 ,1)   npads(1,2, 2 ,2, 4 ,1)   npads(1,2, 3 ,2, 4 ,1);
%           npads(1,2, 1 ,3, 1 ,1)   npads(1,2, 2 ,3, 1 ,1)   npads(1,2, 3 ,3, 1 ,1)   npads(1,2, 1 ,3, 4 ,1)   npads(1,2, 2 ,3, 4 ,1)   npads(1,2, 3 ,3, 4 ,1);
%           npads(1,2, 1 ,4, 1 ,1)   npads(1,2, 2 ,4, 1 ,1)   npads(1,2, 3 ,4, 1 ,1)   npads(1,2, 1 ,4, 4 ,1)   npads(1,2, 2 ,4, 4 ,1)   npads(1,2, 3 ,4, 4 ,1)];
%       
% f1 = figure(1);
% clf
% bar(Yplot,'group')
% set(gca,'yscale','log')
% set(gca,'xtickmode','manual')
% set(gca,'xtick',[1 2 3 4])
% set(gca,'xticklabel',[1 2 4 8])
% xlabel('Number of tiers')
% ylabel('Number of power delivery TSVs')
% fixfigs(1,3,14,12)
% 
% 
% %% Plot the number of PSN TSVs vs. decap
% %npads(cind,dind,thind,nind,pind,freq_ind)
% f2 = figure(2);
% clf
% plot_color = [0 0 0; 0 0 1; 1 0 0; 0 1 0];
% %plot_styles = {'-','-','--','--','-','-','--','--'};
% hold on
% 
% pinds_to_plot = [1 4];
% num_pinds_to_plot = length(pinds_to_plot);
% Yplot2 = [];
% for nind = 1:num_stacks
%     yvec = zeros(1,num_pinds_to_plot*num_decaps);
%     for pindind = 1:num_pinds_to_plot
%         pind = rel_permittivities(pindind);
%         for dind = 1:num_decaps
%             freq_ind = 1;
%             cind = 1;
%             thind = 1;
%             yvec((pindind-1)*num_decaps + dind) = npads(cind,dind,thind,nind,pind,freq_ind);
%         end
%     end
%     
%     %Yplot2(nind,:) = yvec;
%     Yplot2 = [Yplot2; yvec];
% end
%             
% bar(Yplot2,'group')
% set(gca,'yscale','log')
% set(gca,'xtickmode','manual')
% set(gca,'xtick',[1 2 3 4])
% set(gca,'xticklabel',[1 2 4 8])
% xlabel('Number of tiers')
% ylabel('Number of power delivery TSVs')
% fixfigs(2,3,14,12)
% 
% 
% 
% %% Plot the number of PSN TSVs vs. decap
% %npads(cind,dind,thind,nind,pind,freq_ind)
% f3 = figure(3);
% clf
% plot_color = [0 0 0; 0 0 1; 1 0 0; 0 1 0];
% %plot_styles = {'-','-','--','--','-','-','--','--'};
% hold on
% 
% pinds_to_plot = [4];
% num_pinds_to_plot = length(pinds_to_plot);
% Yplot3 = [];
% for nind = 1:num_stacks
%     yvec = zeros(1,num_pinds_to_plot*num_decaps);
%     for pindind = 1:num_pinds_to_plot
%         pind = pinds_to_plot(pindind);
%         for dind = 1:num_decaps
%             freq_ind = 1;
%             cind = 1;
%             thind = 2;
%             yvec((pindind-1)*num_decaps + dind) = npads(cind,dind,thind,nind,pind,freq_ind);
%         end
%     end
%     %yvec
%     %Yplot2(nind,:) = yvec;
%     Yplot3 = [Yplot3; yvec];
% end
%             
% b = bar(Yplot3,0.8,'group');
% set(b(1),'FaceColor',[ 0 0 0.85])
% set(b(2),'FaceColor',[0.9 0.9 0])
% set(b(3),'FaceColor',[0.85 0 0])
% set(gca,'yscale','log')
% set(gca,'xtickmode','manual')
% set(gca,'xtick',[1 2 3 4])
% set(gca,'xticklabel',[1 2 4 8])
% xlabel('Number of tiers')
% ylabel('Number of power delivery TSVs')
% fixfigs(3,3,14,12)
%     







%% Find the maximum frequency that the device can run to stay below 90C
% Tmax = 90;
% fmax_temp_limited = zeros(num_cooling_configs,num_thicks,num_perms,length(tiers));
% 
% for cind = 1:num_cooling_configs
%     for thind = 1:num_thicks
%         for nind = 1:length(tiers)
%             for pind = 1:num_perms
%                 Tind = find( (temp(cind,thind,nind,pind,:) <= Tmax), 1, 'last');
%                 fmax_temp_limited(cind,thind,pind,nind) = frequencies(Tind);
%             end
%         end
%     end
% end


%% Plot frequency chart
% figure(1)
% clf
% %set(gcf,'DefaultAxesColorOrder',[0 0 1; 1 0 0 ])
% %set(gcf,'DefaultAxesLineStyleOrder','-|--|-.|:')
% %hold all
% hold on
% 
% plot_colors = [1 0 0; 0 0 1; 1 0 0; 0 0 1; 1 0 1; 0 1 0; 1 0 1; 0 1 0];
% plot_styles = {'-','-','--','--','-','-','--','--'};
% 
% plot_ind = 1;
% for cind = 1:num_cooling_configs
%     for thind = num_thicks:-1:1
%         for pind = num_perms:-1:1
%             
%             fmax_tl_plot = zeros(1,length(tiers));
%             for nind = 1:length(tiers)
%                 fmax_tl_plot(nind) = fmax_temp_limited(cind,thind,pind,nind);
%             end
%             
%             plot(tiers,fmax_tl_plot,'color',plot_colors(plot_ind,:),'linestyle',plot_styles{plot_ind})
%             plot_ind = plot_ind + 1;
%         end
%     end
% end
% set(gca,'yscale','log')
% %ylim([5e8 0.5e10])
% xlabel('Number of tiers')
% ylabel('Frequency (Hz)')
% fixfigs(1,3,14,12)
% 





