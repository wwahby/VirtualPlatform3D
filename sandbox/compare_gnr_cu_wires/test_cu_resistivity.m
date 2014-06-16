% Test the resistivity functions

wire_thicknesses = (5:500)*1e-9;
grain_sizes = wire_thicknesses;
resistivity_bulk = 17.2e-9;
electron_mfp = 39e-9; % (m) Mean free path of electrons in copper
specularity_coeff = 0.55;
reflection_coeff = 0.43;

rho_cu_vec = zeros(1,length(wire_thicknesses));
delta_rho_fs_vec = zeros(1,length(wire_thicknesses));
delta_rho_ms_vec = zeros(1,length(wire_thicknesses));


for thind = 1:length(wire_thicknesses)
    film_thickness = wire_thicknesses(thind);
    grain_size = grain_sizes(thind);

    [rho_cu delta_rho_fs delta_rho_ms] = cu_resistivity_fsms(resistivity_bulk,film_thickness,grain_size,electron_mfp,specularity_coeff,reflection_coeff);
    rho_cu_vec(thind) = rho_cu;
    delta_rho_fs_vec(thind) = delta_rho_fs;
    delta_rho_ms_vec(thind) = delta_rho_ms;
end

delta_rho_vec = rho_cu_vec - resistivity_bulk;

%% Plots
thick_nm = wire_thicknesses*1e9;

figure(1)
clf
plot(wire_thicknesses*1e9,rho_cu_vec*1e9);
hold on
plot(wire_thicknesses*1e9,resistivity_bulk*1e9*ones(1,length(rho_cu_vec)),'k--')
xlabel('Line width (nm)')
ylabel('Resistivity (Ohm*nm)')
set(gca,'yscale','log')
fixfigs(1,3,14,12)

%%
figure(2)
clf
plot(thick_nm,delta_rho_vec*1e9,'k')
hold on
plot(thick_nm,delta_rho_fs_vec*1e9,'b')
plot(thick_nm,delta_rho_ms_vec*1e9,'r')
xlabel('Line width (nm)')
ylabel('\Delta\rho ({\Omega}nm)')
set(gca,'yscale','log')
fixfigs(2,3,14,12)


%%
figure(4)
clf
plot(wire_thicknesses*1e9,rho_cu_vec/resistivity_bulk);
xlabel('Line width (nm)')
ylabel('\rho_{cu}/\rho_o')
set(gca,'yscale','log')
set(gca,'xscale','log')
fixfigs(4,3,14,12)