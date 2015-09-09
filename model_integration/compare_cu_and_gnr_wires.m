close all
clear all

%% ITRS2010 9.5nm node params

Rsd_w = 120; % Ohm/micron
Cg_w = 0.5e-15; % F/micron
Lg = 8.9e-9;
Wmin_p = 3*Lg;
Wmin_n = Lg;

driver_size_mult = 1;
Wp = driver_size_mult * Wmin_p;
Wn = driver_size_mult * Wmin_n;

Rsd_p = Rsd_w * Wp;
Cg_p = Cg_w * Wp;
Cg_n = Cg_w * Wn;

Cg = Cg_p + Cg_n;

%% Driver/Load Parameters
R_source = Rsd_p; % (Ohms)
C_source = Cg; % (F)
C_load = C_source; % (F)
R_contact_gnr = 4.3e3; % (Ohms)
epsr_dielectric = 1.85; % (-)

%% SB Test case Driver/Load Parameters
% R_source = 8e3;
% C_source = 1e-15 * 1e6 * Wp;
% C_load = C_source;

%% Interconnect parameters
xc_width = 9.5e-9; % (m)
gnr_widths = [2.5e-9, 5e-9, 10e-9, 20e-9];
xc_length_min = 1e-9;
xc_length_max = 1000e-6;
xc_space = xc_width; % (m)
space_vertical = 20e-9;

%% GNR parameters
mfp_eff_gnr = 0.3e-6;
num_layers_gnr = 10;
rho_interlayer_gnr_cm = 0.3;
rho_interlayer_gnr = rho_interlayer_gnr_cm/1e2;
prob_backscattering = 0.25;
Ef_gnr = 0.2;
temp_K = 300;

%% Cu parameters
resistivity_bulk_cu = 17e-9;
mfp_electron_cu = 39e-9;
specularity_coeff = 0.55;
reflection_coeff = 0.45;
horiz_space = xc_space;
R_contact_cu = 100; % random guess for contact resistance of cu to cu
cu_aspect_ratio = 2;
cu_wire_height = cu_aspect_ratio * xc_width;

%%
num_lengths = 1e2;
num_gnr_widths = length(gnr_widths);
xc_length_vec = linspace(xc_length_min, xc_length_max, num_lengths);

gnr_delay_mat= zeros(num_gnr_widths,num_lengths);
gnr_delay_kumar_mat= zeros(num_gnr_widths,num_lengths);
cu_delay_vec = zeros(1,num_lengths);

C_cu_vec = zeros(1,num_lengths);
C_gnr_mat = zeros(num_gnr_widths,num_lengths);

R_gnr_mat = zeros(num_gnr_widths, num_lengths);
R_gnr_kumar_mat = zeros(num_gnr_widths, num_lengths);
R_cu_vec = zeros(1, num_lengths);



for w_ind = 1:num_gnr_widths
    gnr_width = gnr_widths(w_ind);
    for l_ind = 1:num_lengths
        xc_length = xc_length_vec(l_ind);

        [delay_gnr_rc, R_gnr, C_gnr] = xcm.gnr_calc_delay_rc_NS(R_source, C_source, C_load, R_contact_gnr, epsr_dielectric, space_vertical, gnr_width, xc_length, mfp_eff_gnr, num_layers_gnr, rho_interlayer_gnr, Ef_gnr, temp_K);
        [delay_gnr_rc_kumar, R_gnr_kumar, C_gnr_kumar] = xcm.gnr_calc_delay_rc_kumar(R_source, C_source, C_load, R_contact_gnr, epsr_dielectric, space_vertical, gnr_width, xc_length, mfp_eff_gnr, prob_backscattering, num_layers_gnr, rho_interlayer_gnr, Ef_gnr, temp_K);
        [tau_rc_cu, R_wire_cu, C_wire_cu, rho_cu, C_wire_cu_pul] = xcm.calc_cu_wire_rc_const(resistivity_bulk_cu, xc_width, cu_wire_height, xc_length, mfp_electron_cu, specularity_coeff, reflection_coeff, epsr_dielectric, xc_space, space_vertical);
        delay_cu_rc = xcm.calc_delay_elmore( R_source, C_source, C_load, R_contact_cu, R_wire_cu, C_wire_cu);

        gnr_delay_mat(w_ind,l_ind) = delay_gnr_rc;
        gnr_delay_kumar_mat(w_ind,l_ind) = delay_gnr_rc_kumar;
        cu_delay_vec(l_ind) = delay_cu_rc;
        C_cu_vec(l_ind) = C_wire_cu;
        C_gnr_mat(w_ind,l_ind) = C_gnr;
        
        R_gnr_mat(w_ind, l_ind) = R_gnr;
        R_gnr_kumar_mat(w_ind, l_ind) = R_gnr_kumar;
        R_cu_vec(w_ind, l_ind) = R_wire_cu;
    end
end

% Energy consumption
Vdd = 1;
E_cu_vec = 1/2 * (C_source + C_load + C_cu_vec) * Vdd^2;
E_gnr_mat = 1/2* ( C_source + C_load + C_gnr_mat) * Vdd^2;

% EDP
EDP_cu_vec = (E_cu_vec*1e15) .* (cu_delay_vec*1e12); % fJ * ps
EDP_gnr_mat = (E_gnr_mat*1e15) .* (gnr_delay_mat*1e12); % fJ * ps
EDP_gnr_kumar_mat = (E_gnr_mat*1e15) .* (gnr_delay_kumar_mat*1e12); % fJ * ps

%% Figures

colors = [ 0 0 1; 0 1 0; 1 0 0; 1 0 1; 0 1 1; 1 1 0]; 

figure(1)
clf
hold on
plot(xc_length_vec*1e9, cu_delay_vec*1e12, 'k-')
for w_ind = 1:num_gnr_widths
    plot(xc_length_vec*1e9, gnr_delay_mat(w_ind,:)*1e12, 'color', colors(w_ind,:) )
    plot(xc_length_vec*1e9, gnr_delay_kumar_mat(w_ind,:)*1e12, 'color', colors(w_ind,:), 'linestyle', '--' )
end
xlabel('Interconnect Length (nm)')
ylabel('Delay (ps)')
set(gca,'xscale','log')
set(gca,'yscale','log')
fixfigs(1,3,14,12)

figure(2)
clf
hold on
plot(2*xc_length_vec/xc_width, cu_delay_vec*1e12, 'k:')
for w_ind = 1:num_gnr_widths
    plot(2*xc_length_vec/xc_width, gnr_delay_mat(w_ind,:)*1e12, 'color', colors(w_ind,:) )
    plot(2*xc_length_vec/xc_width, gnr_delay_kumar_mat(w_ind,:)*1e12, 'color', colors(w_ind,:), 'linestyle', '--' )
end
xlabel('Interconnect Length (Min Pitches)')
ylabel('Delay (ps)')
set(gca,'xscale','log')
set(gca,'yscale','log')
fixfigs(2,3,14,12)

figure(3)
clf
hold on

for w_ind = 1:num_gnr_widths
    plot(xc_length_vec*1e9, C_gnr_mat(w_ind,:)*1e15, 'color', colors(w_ind,:) )
end
plot(xc_length_vec*1e9, C_cu_vec*1e15,'k:')
xlabel('Interconnect Length (nm)')
ylabel('Capacitance (fF)')
set(gca,'xscale','log')
set(gca,'yscale','log')
fixfigs(3,3,14,12)

figure(4)
clf
hold on

for w_ind = 1:num_gnr_widths
    plot(xc_length_vec*1e9, E_gnr_mat(w_ind,:)*1e15, 'color', colors(w_ind,:) )
end
plot(xc_length_vec*1e9, E_cu_vec*1e15, 'k:')
xlabel('Interconnect Length (nm)')
ylabel('EPB (fJ)')
%set(gca,'xscale','log')
set(gca,'yscale','log')
fixfigs(4,3,14,12)


figure(5)
clf
hold on
plot(xc_length_vec*1e9, EDP_cu_vec, 'k')
for w_ind = 1:num_gnr_widths
    plot(xc_length_vec*1e9, EDP_gnr_mat(w_ind,:), 'color', colors(w_ind,:) )
    plot(xc_length_vec*1e9, EDP_gnr_kumar_mat(w_ind,:), 'color', colors(w_ind,:), 'linestyle', '--' )
end
xlabel('Interconnect Length (nm)')
ylabel('EDP (fJ*ps)')
set(gca,'xscale','log')
set(gca,'yscale','log')
fixfigs(5,3,14,12)

figure(6)
clf
hold on
plot(xc_length_vec*1e9, R_cu_vec/1e3, 'k')
for w_ind = 1:num_gnr_widths
    plot(xc_length_vec*1e9, R_gnr_mat(w_ind,:)/1e3, 'color', colors(w_ind,:) )
    plot(xc_length_vec*1e9, R_gnr_kumar_mat(w_ind,:)/1e3, 'color', colors(w_ind,:), 'linestyle', '--' )
end
xlabel('Interconnect Length (nm)')
ylabel('Resistance (kOhm))')
set(gca,'xscale','log')
set(gca,'yscale','log')
fixfigs(6,3,14,12)



%% 

gnr_width = 8.9e-9;
pitch_min = 2*gnr_width;
xc_length = 100*pitch_min;
R_contact_gnr = 4.3e3;

num_layers_vec = 1:10;
nmax = length(num_layers_vec);

gnr_delay_vec_layers = zeros(1,nmax);
gnr_delay_kumar_vec_layers = zeros(1,nmax);
gnr_delay_kumar_esc_vec_layers = zeros(1,nmax);
cu_delay_vec_layers = zeros(1,nmax);
C_cu_vec_layers = zeros(1,nmax);
C_gnr_vec_layers = zeros(1,nmax);
C_gnr_esc_vec_layers = zeros(1,nmax);
R_gnr_vec_layers = zeros(1,nmax);
R_gnr_kumar_vec_layers = zeros(1,nmax);
R_cu_vec_layers = zeros(1,nmax);

for nind = 1:nmax
    num_layers_gnr = num_layers_vec(nind);
    gnr_height = 0.35e-9 *(2*num_layers_gnr - 1);
    
    [delay_gnr_rc, R_gnr, C_gnr, Rp] = xcm.gnr_calc_delay_rc_NS(R_source, C_source, C_load, R_contact_gnr, epsr_dielectric, space_vertical, gnr_width, xc_length, mfp_eff_gnr, num_layers_gnr, rho_interlayer_gnr, Ef_gnr, temp_K);
    [delay_gnr_rc_kumar, R_gnr_kumar, C_gnr_kumar] = xcm.gnr_calc_delay_rc_kumar(R_source, C_source, C_load, R_contact_gnr, epsr_dielectric, space_vertical, gnr_width, xc_length, mfp_eff_gnr, prob_backscattering, num_layers_gnr, rho_interlayer_gnr, Ef_gnr, temp_K);
    [tau_rc_cu, R_wire_cu, C_wire_cu, rho_cu, C_wire_cu_pul] = xcm.calc_cu_wire_rc_const(resistivity_bulk_cu, xc_width, cu_wire_height, xc_length, mfp_electron_cu, specularity_coeff, reflection_coeff, epsr_dielectric, xc_space, space_vertical);
    delay_cu_rc = xcm.calc_delay_elmore( R_source, C_source, C_load, R_contact_cu, R_wire_cu, C_wire_cu);
    
    % get ESC gnr capacitance
    [C_gnr_esc, C_gnr_esc_pul] = xcm.calc_cu_wire_capacitance(epsr_dielectric, xc_length, gnr_width, gnr_height, xc_space, space_vertical);
    delay_gnr_esc_rc = xcm.calc_delay_elmore( R_source, C_source, C_load, Rp, R_gnr_kumar, C_gnr_esc);
    
    gnr_delay_vec_layers(nind) = delay_gnr_rc;
    gnr_delay_kumar_vec_layers(nind) = delay_gnr_rc_kumar;
    gnr_delay_kumar_esc_vec_layers(nind) = delay_gnr_esc_rc;
    cu_delay_vec_layers(nind) = delay_cu_rc;
    
    C_cu_vec_layers(nind) = C_wire_cu;
    C_gnr_vec_layers(nind) = C_gnr;
    C_gnr_esc_vec_layers(nind) = C_gnr_esc;

    R_gnr_vec_layers(nind) = R_gnr;
    R_gnr_kumar_vec_layers(nind) = R_gnr_kumar;
    R_cu_vec_layers(nind) = R_wire_cu;
end

figure(7)
clf
hold on
plot(num_layers_vec, cu_delay_vec_layers*1e12, 'k:')
plot(num_layers_vec, gnr_delay_kumar_vec_layers*1e12, 'b')
plot(num_layers_vec, gnr_delay_kumar_esc_vec_layers*1e12, 'g')
plot(num_layers_vec, gnr_delay_vec_layers*1e12, 'r-')
xlabel('GNR Layers')
ylabel('Delay (ps)')
fixfigs(7,3,14,12)


%%
delay_target = 0.1e-9;
delay_tolerance = 0.01;
guess_init = 100e-9;
repeater_fraction = 1;
xc_length = 10e-6;
[width_rc, delay_rc, width_rep, delay_rep] = xcm.find_narrowest_gnr_interconnects(delay_target, delay_tolerance, guess_init, repeater_fraction, R_source, C_source, C_load, R_contact_gnr, epsr_dielectric, space_vertical, xc_length, mfp_eff_gnr, num_layers_gnr, rho_interlayer_gnr, Ef_gnr, temp_K);
    
