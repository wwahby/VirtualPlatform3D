%close all

Beta = 0.9;
Ro = 10e3;
Co = 1e-15;
repeater_fraction = [0.5 0.3 0.2];
%repeater_fraction = 0.5;
routing_efficiency_vec = [0.2 0.4];

% [Ln_vec_td pn_vec_td A_wires_td A_vias_wiring_td A_vias_repeaters_td A_layer_td repeater_num_td repeater_size_td tau_rc_vec_td tau_rep_vec_td] = ...
%     wla_topdown_with_repeaters(...
%     iidf_2d,gate_pitch,min_pitch,layers_per_tier,routing_efficiency_vec,...
%     layer_area,rho_m,epsr_d,Beta,Tclk,Rc,Ro,Co,repeater_fraction);

%% pack inputs
chip.iidf = iidf_2d;
chip.gate_pitch = gate_pitch;
chip.min_pitch = min_pitch;
wire.layers_per_tier = layers_per_tier;
wire.routing_efficiency = routing_efficiency_vec;
wire.layer_area = layer_area;
wire.resistivity = rho_m;
wire.dielectric_epsr = epsr_d;
wire.delay_constant = alpha_t;
wire.Beta = Beta;
chip.clock_period = Tclk;
wire.Rc = Rc;
chip.lengths = l2d;
chip.Ro = Ro;
chip.Co = Co;
wire.repeater_fraction = repeater_fraction;
chip.area_total = Ach_m2;

%% Run WLA

wire = wla_topdown_with_repeaters(chip,wire);

%% unpack outputs
Ln_vec_td = wire.Ln;
pn_vec_td = wire.pn;
pn_orig_vec_td = wire.pn_orig;
A_wires_td = wire.wire_area;
A_vias_wiring_td = wire.via_area_wires;
A_vias_repeaters_td = wire.via_area_repeaters;
A_layer_td = wire.area_per_layer;
tau_rc_vec_td = wire.delay_rc;
tau_rep_vec_td = wire.delay_repeaters;
repeater_num_td = wire.repeater_num;
repeater_size_td = wire.repater_size;



%% figures
figure(1)
clf
plot(A_wires_td./A_layer_td,'b')
hold on
plot(A_vias_wiring_td./A_layer_td,'g')
plot(A_vias_repeaters_td./A_layer_td,'r')
fixfigs(1,3,14,12)

figure(2)
clf
plot(pn_vec_td*1e9)
xlabel('metal layer')
ylabel('Wiring pitch (nm)')
fixfigs(2,3,14,12)

figure(3)
clf
plot(repeater_num_td)
xlabel('wire length (GP)')
ylabel('Number of repeaters per interconnect')
fixfigs(3,3,14,12)
