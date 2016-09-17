% {Vp, fft_sys, f, v, Omg, NFFT}
colors = {'b', 'r', 'g', 'k'};
linestyles = {'-', '--'};

thind = 1;
figure(1)
clf
hold on

figure(2)
clf
hold on

figure(3)
clf
hold on

chip_temps = zeros(num_thicks, num_stacks,num_wire_resistivities);
for nind = 1:num_stacks
    for wire_res_ind = 1:num_wire_resistivities
        col_ind = mod(nind-1, length(colors))+1;
        linestyle_ind = mod(wire_res_ind-1, length(linestyles)) + 1;
        psn2 = psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind};
        chip_temps(thind, nind, wire_res_ind) = temperature(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind,num_gates_ind);

        Vp = psn2.output_cell{1};
        fft_sys = psn2.output_cell{2};
        f = psn2.output_cell{3};
        v = psn2.output_cell{4};
        Omg = psn2.output_cell{5};
        NFFT = psn2.output_cell{6};
        Y = psn2.output_cell{7};
        output = psn2.output_cell{8};

        Npts=1024*10; %[FIX] May need a dynamic way to determine FFT resolution. WWAHBY 2014.10.06
        tstop=1000e-9;
        tstep=tstop/Npts;

        figure(1)
        plot(f/1e9,abs(Vp), 'color', colors{col_ind}, 'linestyle', linestyles{linestyle_ind})
        xlabel('Frequency (GHz)')
        ylabel('Voltage Noise')
        set(gca,'yscale','log')
        fixfigs(1,3,14,12)

        opval=ifft(Vp,NFFT,'symmetric');
        max_PSN=abs(min(opval))+ mean(sqrt((opval(1:round(90e-9/tstep))).^2));

        tvec = (0:length(opval)-1)*tstep;
        figure(2)
        plot(tvec(1:100)*1e9,opval(1:100), 'color', colors{col_ind}, 'linestyle', linestyles{linestyle_ind})
        xlabel('Time (ns)')
        ylabel('Voltage Noise')
        fixfigs(2,3,14,12)

        opval=ifft(output,NFFT,'symmetric');
        max_PSN=abs(min(opval))+ mean(sqrt((opval(1:round(90e-9/tstep))).^2));
        tvec = (0:length(opval)-1)*tstep;
        figure(3)
        plot(tvec(1:1200)*1e9,opval(1:1200), 'color', colors{col_ind}, 'linestyle', linestyles{linestyle_ind})
        xlabel('Time (ns)')
        ylabel('Voltage Noise')
        fixfigs(3,3,14,12)
    end
end


thind = 1;
figure(4)
clf
hold on

figure(5)
clf
hold on

figure(6)
clf
hold on
for nind = 1:num_stacks
    for wire_res_ind = 1:num_wire_resistivities
        col_ind = mod(nind-1, length(colors))+1;
        linestyle_ind = mod(wire_res_ind-1, length(linestyles)) + 1;
        psn2 = psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind};
        chip_temps(thind, nind, wire_res_ind) = temperature(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind,num_gates_ind);

        Vp = psn2.output_cell{1};
        fft_sys = psn2.output_cell{2};
        f = psn2.output_cell{3};
        v = psn2.output_cell{4};
        Omg = psn2.output_cell{5};
        NFFT = psn2.output_cell{6};
        Y = psn2.output_cell{7};
        output = psn2.output_cell{8};

        Npts=1024*10; %[FIX] May need a dynamic way to determine FFT resolution. WWAHBY 2014.10.06
        tstop=1000e-9;
        tstep=tstop/Npts;

        figure(4)
        plot(f/1e9,abs(Vp), 'color', colors{col_ind}, 'linestyle', linestyles{linestyle_ind})
        xlabel('Frequency (GHz)')
        ylabel('Voltage Noise')
        set(gca,'yscale','log')
        fixfigs(4,3,14,12)

        opval=ifft(Vp,NFFT,'symmetric');
        max_PSN=abs(min(opval))+ mean(sqrt((opval(1:round(90e-9/tstep))).^2));

        tvec = (0:length(opval)-1)*tstep;
        figure(5)
        plot(tvec(1:100)*1e9,opval(1:100), 'color', colors{col_ind}, 'linestyle', linestyles{linestyle_ind})
        xlabel('Time (ns)')
        ylabel('Voltage Noise')
        fixfigs(5,3,14,12)

        opval=ifft(output,NFFT,'symmetric');
        max_PSN=abs(min(opval))+ mean(sqrt((opval(1:round(90e-9/tstep))).^2));
        tvec = (0:length(opval)-1)*tstep;
        figure(6)
        plot(tvec(1:1200)*1e9,opval(1:1200), 'color', colors{col_ind}, 'linestyle', linestyles{linestyle_ind})
        xlabel('Time (ns)')
        ylabel('Voltage Noise')
        fixfigs(6,3,14,12)
    end
end

thind = 1;
figure(7)
clf
hold on
grid on
for wire_res_ind = 1:num_wire_resistivities
    plot(tiers, squeeze(chip_temps(thind, :, wire_res_ind) ) );
end
ylim([35,70])
fixfigs(7,2,14,12)

thind = 1;
figure(8)
clf
hold on
for wire_res_ind = 1:num_wire_resistivities
    plot(tiers, squeeze(chip_temps(thind, :, wire_res_ind) ) );
end
ylim([35,70])
grid on
fixfigs(8,2,14,12)


figure(9)
clf
hold on
grid on
thind = 1;
wire_res_ind = 1;
plot(tiers, squeeze(chip_temps(thind, :, wire_res_ind)), 'b');

thind = 1;
wire_res_ind = 2;
plot(tiers, squeeze(chip_temps(thind, :, wire_res_ind)), 'r');
set(gca, 'xtick', 1:num_stacks);
set(gca, 'xticklabel', tiers)
fixfigs(9,2,14,12)


psn_pads = zeros(num_thicks, num_stacks, num_wire_resistivities);
tsv_widths = zeros(num_thicks, num_stacks, num_wire_resistivities);
area_per_die_mm2_mat = zeros(num_thicks, num_stacks, num_wire_resistivities);
area_tot_mm2_mat = zeros(num_thicks, num_stacks, num_wire_resistivities);
for thind = 1:num_thicks
    for nind = 1:num_stacks
        for wire_res_ind = 1:num_wire_resistivities
            psn_pads(thind, nind, wire_res_ind) = npads(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind,num_gates_ind);
            tsv_widths(thind, nind, wire_res_ind) = psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind}.power_tsv_width;
            area_per_die_mm2_mat(thind, nind, wire_res_ind) = area_per_layer_mat(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind,num_gates_ind)*1e6;
            area_tot_mm2_mat(thind, nind, wire_res_ind) = area_tot_mat(cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind,num_gates_ind)*1e6;
        
        end
    end
end


% figure(10)
% clf
% hold on
% grid on
% 
% thind = 1;
% wire_res_ind = 2;
% plot(tiers, squeeze(psn_pads(thind, :, wire_res_ind)), 'b');
% 
% thind = 1;
% wire_res_ind = 1;
% plot(tiers, squeeze(psn_pads(thind, :, wire_res_ind)), 'g');
% 
% thind = 3;
% wire_res_ind = 1;
% plot(tiers, squeeze(psn_pads(thind, :, wire_res_ind)), 'r');
% set(gca, 'xtick', 1:num_stacks);
% set(gca, 'xticklabel', tiers)
% fixfigs(10,2,14,12)



tsv_max_area_fraction = 0.01;

psn_tsv_area_vec_m2 = psn_pads.*tsv_widths.^2;
psn_tsv_area_vec_mm2 = psn_tsv_area_vec_m2*1e6;

total_tsv_area_allocation_mm2 = area_per_die_mm2_mat .* tsv_max_area_fraction;
tsv_area_utilization = psn_tsv_area_vec_mm2 ./ total_tsv_area_allocation_mm2;
tsv_area_fraction = psn_tsv_area_vec_mm2 ./ area_per_die_mm2_mat;

% figure(11)
% clf
% hold on
% grid on
% 
% thind = 1;
% wire_res_ind = 2;
% plot(tiers, squeeze(tsv_area_utilization(thind, :, wire_res_ind)), 'b');
% 
% thind = 1;
% wire_res_ind = 1;
% plot(tiers, squeeze(tsv_area_utilization(thind, :, wire_res_ind)), 'g');
% 
% thind = 3;
% wire_res_ind = 1;
% plot(tiers, squeeze(tsv_area_utilization(thind, :, wire_res_ind)), 'r');

set(gca, 'xtick', 1:num_stacks);
set(gca, 'xticklabel', tiers)
set(gca,'yscale','log')
xlabel('Number of tiers')
ylabel('Power TSV Utilization Fraction')
fixfigs(11,2,14,12)

%%
figure(21)
clf
hold on
grid on
linestyles = {'-', '--', ':'};
for nind = 2:num_stacks
    col_ind = mod(nind-1, length(colors));
    for thind = 1:num_thicks
        linestyle_ind = mod(thind-1, length(linestyles)) + 1;
        plot(wire_resistivities*1e9, squeeze(tsv_area_utilization(thind, nind, :)), 'color', colors{col_ind}, 'linestyle', linestyles{linestyle_ind});
    end
end
set(gca,'yscale','log')
xlabel('Wire resistivity (\Omeganm)')
ylabel('Power TSV Utilization Fraction')
fixfigs(21,2,14,12)

%%
figure(12)
clf
hold on
grid on

thind = 1;
wire_res_ind = 2;
plot(tiers, squeeze(tsv_area_fraction(thind, :, wire_res_ind)), 'b');

thind = 1;
wire_res_ind = 1;
plot(tiers, squeeze(tsv_area_fraction(thind, :, wire_res_ind)), 'g');

thind = 3;
wire_res_ind = 1;
plot(tiers, squeeze(tsv_area_fraction(thind, :, wire_res_ind)), 'r');

set(gca, 'xtick', 1:num_stacks);
set(gca, 'xticklabel', tiers)
%xlim([2 4])
set(gca,'yscale','log')
xlabel('Number of tiers')
ylabel('TSV Area Fraction')
fixfigs(12,2,14,12)


barvec = zeros(num_stacks, 3);
thind = 1;
wire_res_ind = 2;
% plot(tiers, squeeze(tsv_area_utilization(thind, :, wire_res_ind)), 'b');
barvec(:, 1) = squeeze(tsv_area_utilization(thind, :, wire_res_ind));

thind = 1;
wire_res_ind = 1;
% plot(tiers, squeeze(tsv_area_utilization(thind, :, wire_res_ind)), 'g');
barvec(:, 2) = squeeze(tsv_area_utilization(thind, :, wire_res_ind));

thind = 3;
wire_res_ind = 1;
% plot(tiers, squeeze(tsv_area_utilization(thind, :, wire_res_ind)), 'r');
barvec(:, 3) = squeeze(tsv_area_utilization(thind, :, wire_res_ind));


%%

figure(13)
clf
hold on
grid on
b = bar(barvec(2:end, :), 1);
set(gca, 'xtick', 1:num_stacks);
set(gca, 'xticklabel', 2:4)
%xlim([2 4])
set(gca,'yscale','log')
xlabel('Number of tiers')
ylabel('Fraction of TSV area used')
fixfigs(13,2,14,12)

b(1).FaceColor = 'b';
b(2).FaceColor = 'g';
b(3).FaceColor = 'y';



Avias = barvec;
Arem = 1 - Avias;
dvias = [0.01 0.1 20];
dvias = [dvias; dvias; dvias; dvias];
Aper = 1/4*pi*(dvias*1e-3).^2;
Apow = Aper.*barvec;
Nv = Arem./Aper;

figure(14)
clf
hold on
grid on
b = bar(Nv(2:end, :), 1);
set(gca, 'xtick', 1:num_stacks);
set(gca, 'xticklabel', 2:4)
%xlim([2 4])
set(gca,'yscale','log')
set(gca, 'ytick', 10.^[2:12])
xlabel('Number of tiers')
ylabel('Signal TSVs remaining')
fixfigs(14,2,14,12)

b(1).FaceColor = 'b';
b(2).FaceColor = 'g';
b(3).FaceColor = 'y';


