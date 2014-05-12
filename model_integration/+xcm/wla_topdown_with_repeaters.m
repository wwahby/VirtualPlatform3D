%function [Ln_vec pn_vec A_wires A_vias_wiring A_vias_repeaters A_layer repeater_num repeater_size tau_rc_vec tau_rep_vec] = wla_topdown_with_repeaters(Iidf,gate_pitch,min_pitch,layers_per_tier,routing_efficiency_vec,layer_area,rho_m,epsr_d,Beta,Tclk,Rc,Ro,Co,repeater_fraction)
function [wire repeater] = wla_topdown_with_repeaters(chip,gate,wire)

% Simple wire layer assignment algorithm, including via blockage
% Assign wires to layers to satisfy timing and areal constraints as
% described by Venkatesan (see Joyner thesis)
% Assuming wires have unity aspect ratio (i,e, they are pn/2 wide and pn/2
% tall)
% via blockage is included by counting how many wires remain to be routed,
% and assuming a square (pn/2)^2 via will be required on wiring layer n for
% each wire not yet routed.
%==
% Inputs
% ==
% Iidf - wirelength distribution (number of wires of a particular length
%        Wire lengths should be in integer units of gate pitches
% gate_pitch - (m) Distance between logic gates (NOT the width of the gate
%              on an individual transistor!)
% min_pitch - (m) minimum patternable feature pitch
% layers_per_tier - (-) how many metal layers per wiring tier? Wires are
%                   often routed in horizontal and vertical tiers, so
%                   Joyner and Venkatesan use 2 for this (since they want
%                   the dimensions of wires in the H and V tiers to be the
%                   same)
% routing_efficiency - (-) Fraction [0-1] of layer area that your routing tool
%                        actually uses for wire routing. Usually ~0.4-0.6
% layer_area - (m^2) Area available in each layer for routing wires
% rho_m - (Ohm*m) Resistivity of wires
%         If a 1xn vector is input, the first n entries will be used for
%         the first n wiring tiers, and all tiers above n will use the nth
%         entry.
% epsr_d - (-) relative permittivity of interlayer dielectric
% alpha_t - (-) timing constant related to wiring capacitance and delay
%            Venkatesan gives this as 1.1*6.2
% Beta - (-) Fraction [0-1] of clock period that wire delay in each
%         tier can consume. Handles vectors the same way rho_m does
% Tclk - (s) Clock period
% Rc - (Ohm) Contact resistance between tiers. Use 0 to ignore Rc
%          Like rho_m, you can input a vector to have different Rc's for
%          different tiers

%% Unpack inputs from input object
Iidf = chip.iidf;
gate_pitch = chip.gate_pitch;
min_pitch = chip.min_pitch;
layers_per_tier = wire.layers_per_tier;
routing_efficiency_vec = wire.routing_efficiency;
layer_area = wire.layer_area;
rho_m = wire.resistivity;
epsr_d = wire.dielectric_epsr;
cap_const = wire.capacitance_constant;
Beta = wire.Beta;
Tclk = chip.clock_period;
Rc = wire.Rc;
Ro = gate.output_resistance;
Co = gate.capacitance;
repeater_fraction = fliplr(wire.repeater_fraction); % flip these around since we start wiring with the top/last metal level in this case


%% Repeater constraints
% Min inverter size (estimate)
Ainv_min = gate_pitch^2/4;
repeater_area_fraction = 0.05;
repeater_via_area_fraction = 0.01;

%%
chi = 2/3; % conversion factor -- converts between point-to-point wirelength and total net length
eps0 = 8.854e-12;
eps_d = epsr_d * eps0;

l = chip.lengths;
Iidf = round(Iidf); % rounding this because it's nonsensical to route 0.32 of a wire
LIDF = l.*Iidf;
lmax = find(Iidf,1,'last'); % find last length for which we have more than zero wires

repeater_num = zeros(1,length(Iidf));
repeater_size = zeros(1,length(Iidf));

% Preallocate a bunch of stuff
% assume we won't use more than some ridiculous number of layers
max_layers = 100;
Ln_vec = zeros(1,max_layers);
pn_vec = zeros(1,max_layers);
pn_orig_vec = zeros(1,max_layers);
A_wires = zeros(1,max_layers);
A_vias_wiring = zeros(1,max_layers);
A_vias_repeaters = zeros(1,max_layers);
A_layer = zeros(1,max_layers);
tau_rc_vec = zeros(1,max_layers);
tau_rc_gnr_vec = zeros(1,max_layers);
tau_rep_vec = zeros(1,max_layers);
material_vec = zeros(1,max_layers);
R_int_vec = zeros(1,max_layers);

n = 1; % start with top wiring tier
Ln = lmax; % length of longest wire routed in this tier (GP)
Lm = lmax + 1; % length of longest wire routed in previous tier (GP) Starts above Lm to ensure we enter the assignment loop below
Ln_vec(1) = lmax;
final_layers = 0;
num_repeater_vias = 0; % No repeater vias through top level
A_vr_n = 0; % No repeater vias through top level
repeater_area_used = 0; %No repeater area used yet
%%

while (Lm >= 0 && n < max_layers)
    Ln_ind = Ln+1;
    Lm_ind = Lm+1;
    
    % Get the value of these constants for each layer we're looking at
    % Each of these can be input as a vector
    % If it's a vector, grab the nth element
    % If n is larger than the last element, grab the last element
    
    rho_m_n = xcm.get_nth_or_last(rho_m,n);
    Rc_n = xcm.get_nth_or_last(Rc,n);
    routing_efficiency = xcm.get_nth_or_last(routing_efficiency_vec,n);
    if( final_layers == 0)
        Beta_n = xcm.get_nth_or_last(Beta,n);
    else
        Beta_n = wire.Beta_short;
    end
    gamma = xcm.get_nth_or_last(repeater_fraction,n);
    
    % Determine layer area available
    A_max_n = layers_per_tier * routing_efficiency * layer_area;
    A_layer(n) = A_max_n;
    
    % Repeater area constraint
    Arep_max = repeater_area_fraction*A_max_n; % [FIX] Need a better way to find available area

    % Repeater via constraint
    Arep_via_max = repeater_via_area_fraction*A_max_n; % [FIX] Need a better way to find available area
    
    Ln_m = @(Ln) gate_pitch*Ln;
    %R_int = @(pn,Ln) 4*rho_m_n*Ln_m(Ln)/pn^2 + Rc_n; % [FIX] This should involve the width fraction and aspect ratio, curretnly assumes wf = 0.5 and ar = 1
    R_int = @(pn,Ln) rho_m_n*Ln_m(Ln)/ (wire.aspect_ratio * wire.width_fraction^2 * pn^2) + Rc_n; % includes impact of differently-sized wires
    C_int = @(Ln) cap_const*eps_d*Ln_m(Ln);
    
    % Fit function to determine delay when using sub-optimal repeater
    % fraction (from Joyner) (gamma)
    alpha_rep = @(gamma) (1.44 + 0.53*(gamma + 1/gamma) );
    
    %pn_rc_min = @(Ln,gamma) sqrt(1.1*cap_const*4*rho_m_n*eps_d/(Beta_n*Tclk - 1.1*cap_const*Rc_n*eps_d*Ln_m(Ln))) * Ln_m(Ln);
    %pn_rep = @(Ln,gamma) max( min_pitch, sqrt(1.1*cap_const*4*rho_m*eps_d / ( (Beta_n*Tclk)^2/(alpha_rep(gamma)^2*Ro*Co) - 1.1*cap_const*Rc*eps_d*Ln_m(Ln) )) * Ln_m(Ln) );
    
    pn_rc_min = @(Ln,gamma) sqrt(1.1*cap_const*rho_m_n*eps_d/(wire.aspect_ratio*wire.width_fraction^2)/(Beta_n*Tclk - 1.1*cap_const*Rc_n*eps_d*Ln_m(Ln))) * Ln_m(Ln);
    pn_rc = @(Ln,gamma) max(min_pitch, pn_rc_min(Ln,gamma) );
    pn_rep = @(Ln,gamma) max( min_pitch, sqrt(1.1*cap_const*4*rho_m*eps_d/(wire.aspect_ratio*wire.width_fraction^2) / ( (Beta_n*Tclk)^2/(alpha_rep(gamma)^2*Ro*Co) - 1.1*cap_const*Rc*eps_d*Ln_m(Ln) )) * Ln_m(Ln) );
    
    % First, figure out what the pitch would be with and without repeaters
    pn_no_rep = pn_rc(Ln,gamma);
    pn_with_rep = pn_rep(Ln,gamma);
    pn_rc_orig = pn_rc_min(Ln,gamma);
    pn_orig_vec(n) = pn_rc_orig;
    R_int_vec(n) = R_int(pn_no_rep,gamma);
    
    % Now figure out what the actual RC time constant of these wires would
    % be. We'll compare this to the min driver delay to determine whether
    % inserting repeaters would be a good idea
    tau_rc = R_int(pn_no_rep,Ln)*C_int(Ln);
    
    % If the RC delay constant is long enough, we'll use repeaters (and use
    % the pitch derived from the repeater delay)
    % Otherwise we'll just use the RC delay to determine the pitch
    use_repeaters = ((tau_rc > 7*Ro*Co) & (A_vr_n < Arep_via_max) & (repeater_area_used < Arep_max) & (pn_with_rep < pn_no_rep ) );
    if(use_repeaters )  % use repeaters
        pn_vec(n) = pn_with_rep;
    else
        pn_vec(n) = pn_no_rep;
    end
    tau_rep = alpha_rep(gamma)*sqrt(Ro*Co*R_int(pn_vec(n),Ln)*C_int(Ln)); % interconnect delay with repeaters
    
    %% Can we do better with graphene?
   
    width_fraction = wire.width_fraction;
    width_cu = pn_vec(n)*width_fraction; % width of actual Cu wires
    % Max allowable delay
    delay_max = Beta_n * Tclk;
    
    % Cu wire delay
    if (use_repeaters)
        delay_cu = tau_rep;
    else
        delay_cu = tau_rc;
    end
    
    % GNR delay
    % [FIX] Should allow user to set these -- add to wire object?
    num_layers = 5;
    gnr_length = Ln_m(Ln); % (m)
    temp_K = chip.temperature+273; % (K) %[FIX] temp not known at this point...
    mfp_defect = 1000e-9; % (m)
    rho_interlayer = 3e-3; % (Ohm cm)
    prob_backscattering = 0; % (-)
    Ef = 0.2; % (eV)
    contact_resistance = 0;
    epsrd = epsr_d;
    height_dielectric = 500e-9; % (m)
    
    % Only consider graphene wires that are smaller than cu wires for now
    % [FIX] Need a better way of picking potential gnr widths
    gnr_widths = [ (4:2:20)*1e-9 (25e-9:10e-9:width_cu) ];
    
    [delay_top_vec delay_side_vec R_top_vec R_top_alt_vec R_side_vec L_vec C_gnr_vec C_gnr_raw_vec Nch_vec mfp_eff_vec ] = ...
        xcm.calc_gnr_params_combined_multiple_widths( ...
        num_layers, gnr_widths, gnr_length, temp_K, mfp_defect, rho_interlayer, prob_backscattering, ...
        Ef, contact_resistance, epsrd, height_dielectric );
    
    % find smallest graphene pitch that beats copper
    last_valid_gnr_ind = find( (delay_top_vec < delay_cu), 1, 'first');
    if (~isempty(last_valid_gnr_ind))
        width_gnr = gnr_widths(last_valid_gnr_ind);
        delay_gnr = delay_top_vec(last_valid_gnr_ind);

        pn_vec(n) = 2*width_gnr; % assuming GNR width fraction is 0.5 for now
        tau_rc_gnr_vec(n) = delay_gnr;

        % Use graphene, no repeaters yet
        material_vec(n) = 2; % GNR
        use_repeaters = 0;
        %tau_rc_gnr =  % Need to update the delay for the wire
    else
        material_vec(n) = 1; % Cu
        tau_rc_gnr_vec(n) = delay_top_vec(end);
    end

    

    %% Now we need to figure out the smallest wire we can route on this tier
    % This is a straightforward application of several factors
    % A_avail_n. Area available for routing
    % Avw_n. Area required for vias to wires on higher levels
    % Avr_n. Area required for repeater vias to higher levels
    % Aw_n. 3. Area required to actually route wires on this layer
    % Need Aw_n + Avr_n + Avw_n <= A_avail_n
    % Avr_n, Avw_n are completely determined by the layers above and the
    % repeater sizing for this layer
    
    via_area_n = (1.5*pn_vec(n))^2; % 1.5 because min via pitch is usually 3*F
 
    % Need to figure out how many vias are required to connect layers above
    % this one and the gates below
    num_repeater_vias = sum( Iidf(Ln+2:end).*repeater_num(Ln+2:end) );
    num_wire_vias = sum( Iidf(Ln+2:end) );
    
    A_vw_n = via_area_n * num_wire_vias;
    A_vr_n = via_area_n * num_repeater_vias;
    
    A_avail_n = A_max_n - A_vw_n - A_vr_n; % This is the area we have remaining for wire routing
    
    % last wire routed in next lowest tier will be the first length for
    % which we don't have enough area
    Areq = 0;
    L_ind = Ln_ind+1;
    Lm_old = Lm;
    while( (Areq < A_avail_n) && (L_ind > 1))
        L_ind = L_ind - 1;
        Areq = Areq + chi*pn_vec(n)*gate_pitch*LIDF(L_ind);
    end
    L = L_ind - 1; % first length routed in this layer
    Lm = L - 1; % last length routed in previous tier

    A_wires(n) = chi*pn_vec(n)*gate_pitch*sum(LIDF(L_ind:Ln_ind));
    A_vias_wiring(n) = A_vw_n;
    A_vias_repeaters(n) = A_vr_n;
    
    %% now that we know how many wires to route in this tier, size repeaters
    wire_lengths_gp = (Lm+1:Ln);
    wire_length_inds = wire_lengths_gp+1;
    
    
    if(use_repeaters) % use repeaters
        % Size repeaters
        repeater_num(wire_length_inds) = gamma*sqrt(0.4*R_int(pn_vec(n),wire_lengths_gp).*C_int(wire_lengths_gp)/0.7/Ro/Co); % number of repeaters
        repeater_size(wire_length_inds) = sqrt(Ro/Co*C_int(wire_lengths_gp)./R_int(pn_vec(n),wire_lengths_gp)); % h*W/L is repeater size
        tau_rep_vec(n) = alpha_rep(gamma)*sqrt(Ro*Co*R_int(pn_vec(n),Ln)*C_int(Ln));
    else
        % Not using repeaters, so just set their number and size to 0
        repeater_num(wire_length_inds) = 0;
        repeater_size(wire_length_inds) = 0;
    end
    tau_rc_vec(n) = tau_rc;
    
    repeater_area_used = 2*Ainv_min*repeater_size(Lm_ind:end).*repeater_num(Lm_ind:end).*Iidf(Lm_ind:end);
    
    if (Lm > 0)
        Ln_vec(n+1) = Lm;
    else
        if(final_layers == 0)
            n = n-1;
            Lm = Lm_old;
            final_layers = 1;
        end
    end
    Ln = Lm;
    
    n = n+1;
end

%% Truncate overallocated vectors
pn_vec = fliplr(pn_vec(Ln_vec > 0));
pn_orig_vec = fliplr(pn_orig_vec(Ln_vec > 0));
Ln_vec = fliplr(Ln_vec(Ln_vec > 0));
A_wires = fliplr(A_wires(Ln_vec > 0));
A_vias_wiring = fliplr(A_vias_wiring(Ln_vec > 0));
A_vias_repeaters = fliplr(A_vias_repeaters(Ln_vec > 0));
A_layer = fliplr(A_layer(Ln_vec > 0));
tau_rc_vec = fliplr(tau_rc_vec(Ln_vec > 0));
tau_rc_gnr_vec = fliplr(tau_rc_gnr_vec(Ln_vec > 0));
tau_rep_vec = fliplr(tau_rep_vec(Ln_vec > 0));
material_vec = fliplr(material_vec(Ln_vec > 0));
R_int_vec = fliplr(R_int_vec(Ln_vec > 0));

%% Pack outputs
wire.Ln = Ln_vec;
wire.pn = pn_vec;
wire.pn_orig = pn_orig_vec;
wire.wire_area = A_wires;
wire.via_area = A_vias_wiring + A_vias_repeaters;
wire.via_area_wires = A_vias_wiring;
wire.via_area_repeaters = A_vias_repeaters;
wire.area_per_layer = A_layer;
wire.delay_rc = tau_rc_vec;
wire.delay_rc_gnr = tau_rc_gnr_vec;
wire.delay_repeaters = tau_rep_vec;
wire.material_vec = material_vec;
wire.resistance_cu = R_int_vec;

[Cxc Cn] = xcm.calc_wiring_capacitance_from_area(wire);
wire.capacitance_total = Cxc;
wire.capacitance_per_tier = Cn;

% clear out any NaNs
repeater_size(isnan(repeater_size)) = 0;
repeater_num(isnan(repeater_num)) = 0;

repeater.num_per_wire = repeater_num;
repeater.size = repeater_size;
repeater.area_total = repeater_area_used;

num_tiers = length(wire.Ln);
repeater.num_per_tier = zeros(1,num_tiers);
tier_start_ind = 1;
for i=1:num_tiers
    Ln_ind = wire.Ln(i)+1;
    repeater.num_per_tier(i) = sum(repeater_num(tier_start_ind:Ln_ind).*Iidf(tier_start_ind:Ln_ind));
    tier_start_ind = Ln_ind+1;
end
    


