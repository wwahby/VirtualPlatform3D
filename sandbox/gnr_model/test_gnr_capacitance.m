%  testing the gnr xc wire capacitance function

%% GNR Capacitance
num_layers = 5;
Ef = 0.2;
temp_K = 300;
widths = (5:100)*1e-9;
gnr_length = 10e-6;
thickness = 0.5e-9*num_layers;
height_dielectric = 300e-9;
epsrd = 3;
eps0 = 8.854e-12;
eps = eps0*epsrd;
prob_backscattering = 0.2;
mfp_defect = 300e-9;

[C_gnr C_gnr_raw cap_const Cqe Calpha] = calc_gnr_wire_capacitance(num_layers,Ef,temp_K,widths,gnr_length,thickness,height_dielectric,epsrd,prob_backscattering,mfp_defect);


%% Cu capacitance
width_fraction = 0.25;
aspect_ratio = 1.8;

wire_length = gnr_length;
pitch = widths/width_fraction;
vertical_spacing = height_dielectric;

[C_cu_scaled C_cu_fixed C_cu_fixed_venk cap_const_cu] = calc_cu_wire_capacitance(epsrd,pitch,wire_length,aspect_ratio,width_fraction,vertical_spacing);
%% Capacitance density per layer
A_wire = (2*widths)*gnr_length; % (m2) 2*w since we're saying the wires are one half-pitch wide
C_m2 = C_gnr_raw./A_wire;
C_mm2 = C_m2/1e6;
C_venk_mm2 = C_gnr./A_wire/1e6;

A_cu_wire = pitch*wire_length; % (m2)
C_cu_m2 = C_cu_fixed./A_cu_wire;
C_cu_mm2 = C_cu_m2/1e6;
C_cu_mm2_venk = C_cu_fixed_venk./A_cu_wire/1e6;

%% Plots
widths_nm = widths*1e9;
figure(1)
clf
plot(widths_nm,C_gnr_raw*1e15,'b')
hold on
plot(widths_nm,C_gnr*1e15,'b--')
plot(widths_nm,C_cu_fixed*1e15,'r')
plot(widths_nm,C_cu_fixed_venk*1e15,'r--')
xlabel('Wire width (nm)')
ylabel('Capacitance (fF)')
set(gca,'yscale','log')
fixfigs(1,3,14,12)

figure(2)
clf
plot(widths_nm,C_mm2)
xlabel('GNR width (nm)')
ylabel('Wiring capacitance density (F/mm^2)')
set(gca,'yscale','log')
fixfigs(2,3,14,12)

figure(3)
clf
plot(widths_nm,C_mm2,'b')
hold on
plot(widths_nm,C_venk_mm2,'b--')
plot(widths_nm,C_cu_mm2,'r')
plot(widths_nm,C_cu_mm2_venk,'r--')
set(gca,'yscale','log')
xlabel('GNR width (nm)')
ylabel('Wiring capacitance density (F/mm^2)')
fixfigs(3,3,14,12)

figure(4)
clf
plot(widths_nm,cap_const)
hold on
plot(widths_nm,cap_const_cu,'r')
xlabel('GNR width (nm)')
ylabel('Venkatesan Number')
fixfigs(4,3,14,12)

figure(5)
clf
plot(widths_nm,cap_const)
hold on
plot(widths_nm,cap_const_cu,'r')
plot(widths_nm,Cqe/eps,'g')
plot(widths_nm,Calpha/eps,'m')
xlabel('GNR width (nm)')
ylabel('C/\epsilon [-]')
set(gca,'yscale','log')
ylim([1e-1 1e3])
fixfigs(5,3,14,12)
