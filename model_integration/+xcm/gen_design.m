%function [iidf l Ln pn pn_orig Cxc Ltot Cn Pdyn Plk Pw Prep Ng_act N_tsvs iidf_rewire] = gen_design(chip,tsv,gate,transistor,wire)
function [chip power wire repeater tsv] = gen_design(chip,tsv,gate,transistor,wire,simulation)

%% Unpack inputs from objects
Ng = chip.num_gates;
alpha = chip.alpha;
k = chip.rent_k;
p = chip.rent_p;
S = chip.num_layers;


Ach_m2 = chip.area_total;
chi = chip.chi;
Tclk = chip.clock_period;
a = chip.logic_activity_factor;
Vdd = chip.Vdd;

gate_pitch = chip.gate_pitch;

w_trans = transistor.gate_length;
eps_ox = transistor.oxide_rel_permittivity;
tox = transistor.oxide_thickness;
Ioff_per_um = transistor.leakage_current_per_micron;

Ro = gate.output_resistance;
N_trans_per_gate = gate.num_transistors;

use_joyner = simulation.use_joyner;
redo_wiring = simulation.redo_wiring_after_repeaters;

AR_tsv = tsv.aspect_ratio;
Atf_max = tsv.max_area_fraction;
h_tsv_m = tsv.height_m;

rho_m = wire.resistivity;
epsr_d = wire.dielectric_epsr;

%% Update some objects where necessary
% [FIX] This assumes the case where all metal layers for the entire 3D
% stack are routed on top of one another (so vias from the top must
% traverse EACH layer)
%wire.layer_area = chip.area_total/chip.num_layers;

% Much more likely is having separate metallization layers near each device
% layer, which significantly reduces the via burden. we can approximate
% this case by just using the entire chip surface area for routing, with
% the understanding that the pitches we find will be defined PER DEVICE
% TIER
% [FIX] Still need a better way to do this. Most likely we'll have to do
% something like this for the first few metal layers, and then when the
% interconnects get long enough we can start considering combined tiers
% Will need to deal with how vias traverse each layer in each case.

if(simulation.separate_wiring_tiers == 1)
    wire.layer_area = chip.area_total; % Determine metal pitch for EACH logic tier
else
    wire.layer_area = Ach_m2/S; % Determine metal pitch for entire wiring stack
end

%% Presize the chip and TSVs
Ns = Ng/S;
Lx = round(sqrt(Ns));
Ach_tier_gp = Ach_m2/gate_pitch^2/S; % Force chip to use specified area 
Ach_tier_m2 = Ach_m2/S;
chip.area_per_layer_gp = Ach_tier_gp;
chip.area_per_layer_m2 = Ach_tier_m2;

% Recalculate these to make sure everything is a nice integer
Ns = Lx^2;
Ng = Ns*S;

% Size the TSVs
h_tsv = ceil(h_tsv_m/gate_pitch);
w_tsv = ceil(h_tsv/AR_tsv);

%% Size the chip so we have a nicely divisible number of unit cells per side
if ((S == 1) || (use_joyner == 1)) % no TSVs for single layer device, or when using original Joyner model
    w_tsv = 0;
    Lxc = Lx;
    Nsc = Ns;
    Ngc = Ng;
    Nuc_1d = 1;
    N_tsvs = 0;
    Tc = 0;
else
    %Nsp = floor( (1+Atf_max)*Ns );
    Nsp = floor(Ns/(1-Atf_max));
    Lxp = floor(sqrt(Nsp));
    Nsp = Lxp^2;
    Tp = ceil(w_tsv/sqrt(Atf_max));

    slack = 0.2;
    [Lxc Tc Nuc_1d gfrac_L gfrac_T] = xcm.find_LT_combination(Lxp,Tp,slack);
    N_tsvs = Nuc_1d^2;
end

Nsc = Lxc^2;
Ngc = Nsc*S;

g_tsv = (Nuc_1d*w_tsv)^2; % number of gates displaced by TSVs
Atf_act = g_tsv/Nsc;

tsv.pitch_gp = Tc;
tsv.pitch_m = Tc*gate_pitch;
tsv.actual_area_fraction = Atf_act;
tsv.num = N_tsvs;

Ns_act = Nsc - g_tsv;
Ng_act = Ns_act*S;
chip.Ng_actual = Ng_act;

% repstr1 = sprintf('Ng_nom: %.4g\tNg_cor: %.4g\tNg_act: %.4g\tAtf_act: %.4g',Ng,Ngc,Ng_act,Atf_act);
% disp(repstr1)

%% Calculate WLD
iidf = xcm.calc_Iidf_corrected(alpha,k,p,Lx,S,h_tsv,Nuc_1d,w_tsv);
%iidf = calc_Iidf(alpha,k,p,round(sqrt(Ng)),1,h_tsv);

%% Cleanup - Get rid of NaNs
iidf(isnan(iidf)) = 0;
lmax = length(iidf) - 1;
l = 0:lmax;

chip.iidf = iidf;
chip.lengths = l;

%% Wire Layer Assignment (WLA) and Repeater Insertion (RI)

wire.capacitance_constant = xcm.calc_capacitance_constant(wire.aspect_ratio,wire.width_fraction);
if (simulation.topdown_WLARI == 1)
    
    % Simultaneously do WLA and RI, starting with TOP metal layer and working
    % downward.
    % This method has the advantage of accurately determining wiring pitch
    % based on the POST-RI wire delay, as well as knowing in advance the total
    % repeater and wire via area required on lower levels. The disadvantage is
    % that the bottom metal layer may not be well-utilized.
    % To fix issues with poor bottom layer fills, we can rerun the routine
    % with different top-level metal fill factors until we find a good
    % wiring solution
    
    [wire repeater] = xcm.wlatdri(chip,gate,wire,simulation);
    
else
    
    % Do sequential WLA and RI. In this method the repeater via area is not
    % known, and the real wire delay (after repeaters) is not known during
    % the initial WLA.
    % These problems can be circumvented by doing this step several times
    % to converge on the real result, but for now we're just doing a
    % single-shot.
    
    % Determine wire pitch and layer assignment
    wire = xcm.wla_improved(chip,wire);

    % Repeater insertion
    repeater = xcm.repeater_insertion(chip,gate,transistor,wire);
end


%% Power estimates
eps0 = 8.854e-12; % (F/m)

thermal_voltage = 0.0258 * (chip.temperature+273)/300; % kT/q (V)
swing_at_temp = transistor.subthreshold_swing * (chip.temperature+273)/300; % (V) Subthreshold swing at the actual chip temperature
Vth = transistor.Vt;
Vgs = 0;
Vds = chip.Vdd;
Vdd = chip.Vdd;
Ilk = (Ioff_per_um*1e6)*(1.5*w_trans);% (1.5 because pmos should have ~3X nmos width, and half the transistors should be pmos)
%[FIX] very coarse leakage model -- update this with something better (incl temp, gate size, etc)

width_um = 1.5*w_trans*1e6;
w_trans_um = w_trans * 1e6;
Ilk_To = Ioff_per_um*width_um;
To = transistor.leakage_reference_temperature;

kb = 1.381e-23; % (J/K)
T = chip.temperature + 273.15; % (K)
q = 1.602e-19; % (C)
Eg = 1.12; % (eV)

phi_th = kb*T/q;
Ilk_T = Ilk_To * (T/To)*exp(-q*Eg/kb*(1/T-1/To))*(exp(-Vds/phi_th)-1);
Ilk_T = abs(Ilk_T);

Plk_per_transistor = Ilk_T*Vdd;
Plk_per_transistor_per_um = Ilk_T/width_um;



%Ct = transistor.capacitance_per_micron * w_trans/1e-6;
%Cox = 4*transistor.capacitance; % (4 because pmos + nmos capacitance in parallel, pmos cap ~3X nmos)
Cox = 1.5*transistor.capacitance; % (1.5 because pmos should have ~3X nmos width, and half the transistors should be pmos)
Co = Cox; % Need to include parasitics for realistic estimate
Co_rep = Cox*repeater.size;
Ilk_rep = Ilk*repeater.size;

Nt = N_trans_per_gate * Ng;
f = 1/Tclk;
Cxc = wire.capacitance_total;

Pdyn = 1/2*a*Co*Vdd^2*f*Nt;
Plk = Vdd*Ilk*Nt;
Plk_logic = Plk_per_transistor * Nt;

Pw = 1/2*a*Cxc*Vdd^2*f;

%Plk_rep_vec = Vdd*Ilk_rep.* repeater.num_per_wire .* repeater.size;
Plk_rep_vec = Plk_per_transistor_per_um .* w_trans_um .* repeater.size .* repeater.num_per_wire .* chip.iidf;
Plk_rep = sum(Plk_rep_vec);

Pdyn_rep_vec = 1/2*a*Vdd^2*f.* repeater.num_per_wire .*Co_rep.* repeater.size;
Pdyn_rep = sum(Pdyn_rep_vec);

Prep = Pdyn_rep + Plk_rep;

Arep_used_mm2 = repeater.area_total*(1e3)^2;

%% update output objects
chip.iidf = iidf;
chip.lengths = l;

power.dynamic = Pdyn;
power.leakage = Plk_logic;
power.wiring = Pw;
power.repeater = Prep;
power.repeater_leakage = Plk_rep;
power.repeater_dynamic = Pdyn_rep;
power.total = Pdyn + Plk_logic + Pw + Prep;
power.density = power.total/chip.area_per_layer_m2;


