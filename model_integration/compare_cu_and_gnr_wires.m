

%% Driver/Load Parameters
R_source = 10e3; % (Ohms)
C_source = 10e-15; % (F)
C_load = C_source; % (F)
R_contact_gnr = 1e3; % (Ohms)
epsr_dielectric = 3; % (-)
%space_vertical = 50e-9;

%% Interconnect parameters
xc_width = 5e-9; % (m)
gnr_widths = [5e-9, 10e-9, 20e-9, 100e-9];
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
num_gnr_widths = length(gnr_widths);
xc_length_vec = linspace(1e-9, 1e-6, num_lengths);
gnr_delay_mat= zeros(num_gnr_widths,num_lengths);
cu_delay_vec = zeros(1,num_lengths);
C_cu_vec = zeros(1,num_lengths);
C_gnr_mat = zeros(num_gnr_widths,num_lengths);

space_vertical = cu_wire_height;

for w_ind = 1:num_gnr_widths
    gnr_width = gnr_widths(w_ind);
    for l_ind = 1:num_lengths
        xc_length = xc_length_vec(l_ind);

        [delay_gnr_rc, R_gnr, C_gnr] = xcm.gnr_calc_delay_rc_NS(R_source, C_source, C_load, R_contact_gnr, epsr_dielectric, space_vertical, gnr_width, xc_length, mfp_eff_gnr, num_layers_gnr, rho_interlayer_gnr, Ef_gnr, temp_K);
        [tau_rc_cu, R_wire_cu, C_wire_cu, rho_cu, C_wire_cu_pul] = xcm.calc_cu_wire_rc_const(resistivity_bulk_cu, xc_width, cu_wire_height, xc_length, mfp_electron_cu, specularity_coeff, reflection_coeff, epsr_dielectric, xc_space, space_vertical);
        delay_cu_rc = xcm.calc_delay_elmore( R_source, C_source, C_load, R_contact_cu, R_wire_cu, C_wire_cu);

        gnr_delay_mat(w_ind,l_ind) = delay_gnr_rc;
        cu_delay_vec(l_ind) = delay_cu_rc;
        C_cu_vec(l_ind) = C_wire_cu;
        C_gnr_mat(w_ind,l_ind) = C_gnr;
    end
end

%% Figures

colors = [ 0 0 1; 0 1 0; 1 0 0; 1 0 1]; 

figure(1)
clf
hold on
plot(xc_length_vec*1e9, cu_delay_vec*1e12, 'k-')
for w_ind = 1:num_gnr_widths
    plot(xc_length_vec*1e9, gnr_delay_mat(w_ind,:)*1e12, 'color', colors(w_ind,:) )
end
xlabel('Interconnect Length (nm)')
ylabel('Delay (ps)')
set(gca,'xscale','log')
set(gca,'yscale','log')
fixfigs(1,3,14,12)

figure(2)
clf
hold on
plot(2*xc_length_vec/xc_width, cu_delay_vec*1e12, 'k-')
for w_ind = 1:num_gnr_widths
    plot(2*xc_length_vec/xc_width, gnr_delay_mat(w_ind,:)*1e12, 'color', colors(w_ind,:) )
end
xlabel('Interconnect Length (Min Pitches)')
ylabel('Delay (ps)')
set(gca,'xscale','log')
set(gca,'yscale','log')
fixfigs(2,3,14,12)

figure(3)
clf
hold on
plot(xc_length_vec*1e9, C_cu_vec*1e15,'k')
for w_ind = 1:num_gnr_widths
    plot(xc_length_vec*1e9, C_gnr_mat(w_ind,:)*1e12, 'color', colors(w_ind,:) )
end
xlabel('Interconnect Length (nm)')
ylabel('Capacitance (fF)')
set(gca,'xscale','log')
set(gca,'yscale','log')
fixfigs(3,3,14,12)