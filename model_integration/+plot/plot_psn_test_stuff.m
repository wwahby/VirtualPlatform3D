% {Vp, fft_sys, f, v, Omg, NFFT}
colors = {'b', 'r', 'g', 'k'};
linestyles = {'-', '--'};

figure(1)
clf
hold on

figure(2)
clf
hold on

figure(3)
clf
hold on
for wire_res_ind = 1:num_wire_resistivities
    col_ind = mod(wire_res_ind-1, length(colors))+1;
    linestyle_ind = mod(wire_res_ind-1, length(linestyles)) + 1;
    psn2 = psn_cell{cind,dind,thind,nind,pind,freq_ind,wire_res_ind,wire_flag_ind,scaling_ind,bar_thick_ind,bar_res_ind,forced_power_ind,k_ind};

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
    plot(f/1e9,abs(Vp))
    xlabel('Frequency (GHz)')
    ylabel('Voltage Noise')
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
    max_PSN=abs(min(opval))+ mean(sqrt((opval(1:round(90e-9/tstep))).^2))
    tvec = (0:length(opval)-1)*tstep;
    figure(3)
    plot(tvec(1:1200)*1e9,opval(1:1200), 'color', colors{col_ind}, 'linestyle', linestyles{linestyle_ind})
    xlabel('Time (ns)')
    ylabel('Voltage Noise')
    fixfigs(3,3,14,12)
end
