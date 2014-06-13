function [wire repeater] = wlatdri(chip,gate,wire)


bottom_layer_underfilled = 1;
min_fill_factor = 0.8; % minimum utilization of available bottom layer area
top_fill_factor = wire.routing_efficiency(1);
wire.routing_efficiency = [wire.routing_efficiency(1) wire.routing_efficiency]; % create a separate entry for first layer

% dead man counter to get out of while loop if something goes wrong
wla_attempts = 0;
max_wla_attempts = 20;
min_top_fill_factor = 0.10;

% Keep track of the best solutions so far
best_bot_fill_factor = 0;
best_top_fill_factor = top_fill_factor;
best_num_wiring_tiers = 9e12; % Absurdly large number so we immediately overwrite this;

% automatically decrease top layer utilization to keep bottom layer
% reasonably filled. Need to do this since TDWLARI doesn't guarantee
% good usage of the bottom layer. Since this routine runs very quickly,
% there isn't really any significant performance overhead, and this
% lets us do a proper job of wire layer assignment and repeater
% insertion while accurately considering wire and repeater via blockage
while((bottom_layer_underfilled == 1) && (wla_attempts < max_wla_attempts) && (top_fill_factor > min_top_fill_factor))
    
    % Actually run topdown WLA and RI
    %[wire_temp repeater_temp] = xcm.wla_topdown_with_repeaters(chip,gate,wire);
    [wire_temp repeater_temp] = xcm.wlatdri_cu_gnr(chip,gate,wire);

    % Figure out what the bottom-layer utilization is
    bot_fill_factor = wire_temp.wire_area(1)/(wire_temp.area_per_layer(1));

    % Figure out whether we need to decrease top layer use, or stop
    if(bot_fill_factor < min_fill_factor)
        top_fill_factor = 0.90*top_fill_factor;
        wire.routing_efficiency(1) = top_fill_factor;
    else
        bottom_layer_underfilled = 0;
    end

    % If we have a better fill factor on the bottom and a smaller
    % number of wiring tiers, this is the best solution
    
    better_wiring_tiers = length(wire_temp.pn) < best_num_wiring_tiers;
    better_bot_fill = bot_fill_factor > best_bot_fill_factor;
    
    
    if( (bot_fill_factor >= best_bot_fill_factor) && (length(wire_temp.pn) <= best_num_wiring_tiers) )
        best_bot_fill_factor = bot_fill_factor;
        best_top_fill_factor = top_fill_factor;
        best_num_wiring_tiers = length(wire_temp.pn);
        wire_best = wire_temp;
        repeater_best = repeater_temp;
        
        fprintf('WLA_BEST: Bot_FF: %.3g \t Top_FF %.3g \t Num_tiers: %.3g\n\n',best_bot_fill_factor,best_top_fill_factor,best_num_wiring_tiers')
    end

    fprintf('WLA: M0 Underfilled: %d \t Top fill factor: %.3g \t Bottom fill factor: %.3g \t Num tiers: %d\n\n',bottom_layer_underfilled,top_fill_factor,bot_fill_factor, length(wire_temp.pn));
    wla_attempts = wla_attempts + 1;
end

% If we still haven't met our target, switch to the best one we found
if(bottom_layer_underfilled == 1)
    fprintf('WLA Warning: Poor M0 fill factor could not be fixed. Reverting to best settings...\n\n')
end
% Revert settings to best found, if we aren't already there
if(best_top_fill_factor ~= top_fill_factor)
    fprintf('WLA: Using best settings...\n\n')
    fprintf('WLA_BEST: Bot_FF: %.3g \t Top_FF %.3g \t Num_tiers: %.3g\n\n',best_bot_fill_factor,best_top_fill_factor,best_num_wiring_tiers')
    wire.routing_efficiency(1) = best_top_fill_factor; % Revert to best fill factor we found previously
    %[wire_temp repeater_temp] = xcm.wlatdri_cu_gnr(chip,gate,wire);
    wire_temp = wire_best;
    repeater_temp = repeater_best;
end

wire = wire_temp;
repeater = repeater_temp;
wire.wla_attempts = wla_attempts;