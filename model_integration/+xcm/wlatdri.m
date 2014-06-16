function [wire repeater] = wlatdri(chip,gate,wire)

bottom_layer_underfilled = 1;
min_fill_factor = 0.9; % minimum utilization of available bottom layer area
top_fill_factor = wire.routing_efficiency(1);
if(length(wire.routing_efficiency) == 1) % if we only have one entry for routing efficiency we need to preserve the nominal value for all the underlying layers
    wire.routing_efficiency = [wire.routing_efficiency(1) wire.routing_efficiency]; % create a separate entry for first layer
end

% dead man counter to get out of while loop if something goes wrong
max_wla_attempts = 16;
min_top_fill_factor = 0.01;

lbnd = min_top_fill_factor;
rbnd = top_fill_factor;

wire.wla_attempts = 2; % 1 because one of the runs below will get thrown away.

wire.routing_efficiency(1) = lbnd;
[wire_lbnd repeater_lbnd] = xcm.wlatdri_cu_gnr(chip,gate,wire);

wire.routing_efficiency(1) = rbnd;
[wire_rbnd repeater_rbnd] = xcm.wlatdri_cu_gnr(chip,gate,wire);


bin_gen = 1;
while( (bottom_layer_underfilled) && (bin_gen <= max_wla_attempts) )
   
    mid = (lbnd + rbnd)/2;
    wire.routing_efficiency(1) = mid;
    
    [wire_mid repeater_mid] = xcm.wlatdri_cu_gnr(chip,gate,wire);
    
    ntiers_lbnd = length(wire_lbnd.pn);
    ntiers_mid = length(wire_mid.pn);
    ntiers_rbnd = length(wire_rbnd.pn);
    
    bot_ff_lbnd = wire_lbnd.wire_area(1)/wire_lbnd.area_per_layer(1);
    bot_ff_mid = wire_mid.wire_area(1)/wire_mid.area_per_layer(1);
    bot_ff_rbnd = wire_rbnd.wire_area(1)/wire_rbnd.area_per_layer(1);
    
    tier_vec = [ntiers_lbnd ntiers_mid ntiers_rbnd];
    bot_ff_vec = [bot_ff_lbnd bot_ff_mid bot_ff_rbnd];
    [min_tiers best_choice_tiers] = min(tier_vec);
    %[min_ff best_choice_ff] = min(bot_ff_vec);
    
    % First, are there multiple choices with the best number of tiers?
    min_tiers_matches = (tier_vec == min_tiers);
    
    bot_ff_vec(~min_tiers_matches) = -1e9; % kill off any choices with higher tier numbers
    
    % Now find best ff from among min_tiers possibilities
    [max_ff best_choice] = max(bot_ff_vec);
    
    if(best_choice == 1) % lbnd
        lbnd = lbnd;
        rbnd = mid;
        
        wire = wire_lbnd;
        repeater = repeater_lbnd;
        
        wire_rbnd = wire_mid;
        repeater_rbnd = repeater_mid;
    elseif(best_choice == 3) %rbnd
        lbnd = mid;
        rbnd = rbnd;
        
        wire = wire_rbnd;
        repeater = repeater_rbnd;
        
        wire_lbnd = wire_mid;
        repeater_lbnd = repeater_mid;
    else % mid is best
        lbnd = lbnd;
        rbnd = mid;
        
        wire = wire_mid;
        repeater = repeater_mid;
        
        wire_rbnd = wire_mid;
        repeater_rbnd = repeater_mid;
    end
    
    best_top_ff = wire.routing_efficiency(1);
    bottom_layer_underfilled = (max_ff < min_fill_factor);
    
    if(bottom_layer_underfilled)
        fprintf('WLA: M0 underfilled! Adjusting top layer fill factor and rerunning...\n')
    end
        
    
    %fprintf('WLA: M0 Underfilled: %d \t Top fill factor: %.3g \t Bottom fill factor: %.3g \t Num tiers: %d \t bin_gen %d\n\n',bottom_layer_underfilled,best_top_ff,max_ff, min_tiers, bin_gen);
    bin_gen = bin_gen + 1;
    
end