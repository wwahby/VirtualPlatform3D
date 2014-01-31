function [Ctot Cn] = calc_wiring_capacitance_from_area_old(pn,layers_per_tier,Awn,Avn,epsrd)
% Calculates total capacitance and capacitance per tier for a wiring
% network, using wiring area per tier as input
% pn - wiring pitch (m) - 1xn vector
% layers_per_tier - how many metal layers are in each wiring tier
% Awn - wiring area on tier n (m^2) Calculated as pn*Lwires_n, so this
%                               includes whitespace because w_wire = pn/2
% Avn - via area on tier n (m^2) Calculated same way as Awn
% epsrd - relative dielectric permittivity of interlayer dielectric

% Calculate permittivity of interlayer dielectric
eps0 = 8.854e-12; % (F/m)
epsd = epsrd*eps0;

num_tiers = length(pn);
Cn = zeros(1,num_tiers);

for i = 1:num_tiers
    
    % dividing Awn by layers per tier to get wiring area per layer within this tier
    % via area does not get divided, since vias go through all layers in each tier
    A_wires_n = Awn(i)/layers_per_tier; % Area devoted to wiring in each layer
    A_vias_n = Avn(i); % Area devoted to vias in each layer
    
    %intra-layer capacitance
    % capacitance between wires in the same layer
    C_mid_wires = epsd/pn(i)*A_wires_n; 
    C_mid_vias = epsd/pn(i)*A_vias_n; 
    C_mid = 2*layers_per_tier* (C_mid_wires + C_mid_vias); % multiplying by 2 because wires have capacitances to either side within a given layer
        
    % no via capacitance in top/bot capacitance calc since vias are just
    % vertical wires -- their capacitance is captured in the intra-layer
    % (C_mid) term
    C_top_a = 1/2*epsd*A_wires_n/pn(i); % capacitance to upper layer for layers below top layer in this wiring tier
    if(i < num_tiers)
        C_top_b = 1/2*epsd*A_wires_n/pn(i+1); % capacitance to upper layer for the top layer in this wiring tier
    else
        C_top_b = 0;
    end
    C_top = (layers_per_tier - 1)*C_top_a + C_top_b;
    
    C_bot_a = 1/2*epsd*A_wires_n/pn(i); % capacitance to bottom layer for layers above bottom layer in this wiring tier
    if(i > 1)
        C_bot_b = 1/2*epsd*A_wires_n/pn(i-1);
    else
        C_bot_b = 0;
    end
    C_bot = (layers_per_tier-1)*C_bot_a + C_bot_b;
    
    C_tier = C_mid + C_bot + C_top; % 2*Cm since we have capacitance on both sides of each wire in the tier
    
    Cn(i) = C_tier;
end

Ctot = sum(Cn);
        
    
     
    
        