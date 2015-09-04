

%% Driver/Load Parameters
R_source = 10e3; % (Ohms)
C_source = 10e-15; % (F)
C_load = C_source; % (F)
R_contact_gnr = 1e3; % (Ohms)
epsr_dielectric = 3; % (-)
%space_vertical = 50e-9;

%% Interconnect parameters
xc_width = 5e-9; % (m)
xc_length = 10e-6; % (m)
xc_space = xc_width; % (m)

%% GNR parameters
mfp_eff_gnr = 1e-6;
num_layers_gnr = 10;
rho_interlayer_gnr_cm = 0.3;
rho_interlayer_gnr = rho_interlayer_gnr_cm/1e2;
Ef_gnr = 0.3;
temp_K = 300;

%% Cu parameters
resistivity_bulk_cu = 17e-9;
mfp_electron_cu = 29e-9;
specularity_coeff = 0.55;
reflection_coeff = 0.45;
horiz_space = xc_space;
R_contact_cu = 100; % random guess for contact resistance of cu to cu
cu_aspect_ratio = 2;
cu_wire_height = cu_aspect_ratio * xc_width;


num_lengths = 1e3;
xc_length_vec = linspace(1e-9, 100e-6, num_lengths);
gnr_delay_vec = zeros(1,num_lengths);
cu_delay_vec = zeros(1,num_lengths);

space_vertical = cu_wire_height;

for l_ind = 1:num_lengths
    xc_length = xc_length_vec(l_ind);

    [delay_gnr_rc, R_gnr, C_gnr] = xcm.gnr_calc_delay_rc_NS(R_source, C_source, C_load, R_contact_gnr, epsr_dielectric, space_vertical, xc_width, xc_length, mfp_eff_gnr, num_layers_gnr, rho_interlayer_gnr, Ef_gnr, temp_K);
    [tau_rc_cu, R_wire_cu, C_wire_cu, rho_cu, C_wire_cu_pul] = xcm.calc_cu_wire_rc_const(resistivity_bulk_cu, xc_width, cu_wire_height, xc_length, mfp_electron_cu, specularity_coeff, reflection_coeff, epsr_dielectric, xc_space, space_vertical);
    delay_cu_rc = xcm.calc_delay_elmore( R_source, C_source, C_load, R_contact_cu, R_wire_cu, C_wire_cu);

    gnr_delay_vec(l_ind) = delay_gnr_rc;
    cu_delay_vec(l_ind) = delay_cu_rc;
end

%% Figures
figure(1)
clf
hold on
plot(xc_length_vec*1e6, cu_delay_vec*1e12, 'k-')
plot(xc_length_vec*1e6, gnr_delay_vec*1e12, 'r')
xlabel('Interconnect Length (microns)')
ylabel('Delay (ps)')
set(gca,'xscale','log')
set(gca,'yscale','log')
fixfigs(1,3,14,12)

figure(2)
clf
hold on
plot(2*xc_length_vec/xc_width, cu_delay_vec*1e12, 'k-')
plot(2*xc_length_vec/xc_width, gnr_delay_vec*1e12, 'r')
xlabel('Interconnect Length (Min Pitches)')
ylabel('Delay (ps)')
set(gca,'xscale','log')
set(gca,'yscale','log')
fixfigs(2,3,14,12)