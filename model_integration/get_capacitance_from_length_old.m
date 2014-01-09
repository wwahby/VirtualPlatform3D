function C_l = get_capacitance_from_length_old(l,Ln_vec,pn_vec,epsr_d,gate_pitch)
% gets capacitance for a wire of length l
% l (gate pitches) length of wire to find capacitance for
% Ln_vec (GP) 1xn vector of longest wire routed in each wiring tier
% pn_vec (m) 1xn vector of wiring pitch on each wiring tier
% epsr_d (-) dielectric constant of interlayer dielectric
% gate_pitch (m) how far gates are spaced from one another

% get dielectric permittivity
eps0 = 8.854e-12; % (F/m) vacuum permittivity
eps = epsr_d*eps0;

wiring_tier = find( l<=Ln_vec,1,'last');
num_tiers = length(Ln_vec);

wire_pitch = pn_vec(wiring_tier);

l_m = l*gate_pitch;

Cm = eps*l_m; % intra-tier capacitance reduces to depend only on wirelength

if (wiring_tier < num_tiers)
    Ct = eps*(l_m/2)*(wire_pitch)/pn(wiring_tier+1);
else
    Ct = eps*(l_m/2); % since we don't know which layer in this tier the wire is routed in, assume worst case (in this case, it's on the bottom of the top tier)
end

Cb = eps*l_m/2; % assume worst case -- wire is in top layer of wiring tier

C_l = 2*Cm + Cb + Ct;
