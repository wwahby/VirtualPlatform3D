%function [Iidf_rep h_vec k_vec Arep_used num_vec size_vec] = repeater_insertion(Iidf,Ach,Ainv_min,pn,Ln,Cn,rho_xcn,Ro_n,Co,w_gate)
function repeater = repeater_insertion(chip,gate,transistor,wire)
%Ach (m^2)
% Ainv_min (m^2)

%% unpack inputs
Iidf = chip.iidf;
Ach = chip.area_total/chip.num_layers;
Ainv_min = 9*chip.gate_pitch^2;
pn = wire.pn;
Ln = wire.Ln;
%Cn = wire.Cn;
rho_xcn = wire.resistivity;
Ro_n = gate.output_resistance;
Co = gate.capacitance;
w_gate = transistor.gate_length;
gate_pitch = chip.gate_pitch;


Iidf = round(Iidf); % let's just deal with integer numbers of interconnects
Iidf = Iidf(2:end); % cut off the zero-length entry
lmax = length(Iidf);

% Find area available for repeaters
%Arep_tot = calc_repeater_area(); % Total area available for repeaters
Arep_tot = 0.10*Ach; % [FIX] Need a better way to find available area
Arep_rem = Arep_tot; % Area remaining for repeaters
Arep_used = 0;

% Conditions for continued repeater insertion
% 1. We still have area available for more repeaters (Arep_rem > 0)
% 2. We still have interconnects of length > 0 (lmax_cur > 0)
% 3. We still get some benefit from inserting repeaters (RxcCxc >= 7RoCo)

add_repeaters = 1;
lmax_cur = find(Iidf>0,1,'last'); % find last nonzero XC
lmax_last = lmax_cur;

h_vec = zeros(1,lmax); % repeater area (compared to min inv)
k_vec = zeros(1,lmax); % number of repeaters per XC
num_vec = zeros(1,length(Ln));
size_vec = zeros(1,length(Ln));

while (add_repeaters == 1)
    
    
    % find properties of the interconnects we're looking at
    % have to multiply lmax_cur by w_gate to get real length
    xc_tier = find(lmax_cur <= Ln,1,'first');
    rho_xc = xcm.get_nth_or_last(rho_xcn,xc_tier);
    Ro = xcm.get_nth_or_last(Ro_n,xc_tier);
    
    %Cxc = Cn(xc_tier)*lmax_cur*w_gate;
    Cxc = xcm.get_capacitance_from_length(lmax_cur,chip,wire); %
    Rxc = rho_xc*lmax_cur*gate_pitch/pn(xc_tier)^2; % pn has units [m]
    
    if(Rxc*Cxc >= 7*Ro*Co)
    % Make sure we actually get a benefit from inserting repeaters
    
        % Size and number the repeaters for each XC
        k = round(sqrt(0.4*Rxc*Cxc/0.7/Ro/Co)); % number of repeaters per XC
        h = sqrt(Ro*Cxc/Rxc/Co); % size of repeaters (compared to min inv)

        num_ins = k*Iidf(lmax_cur);
        Arep_ins = Ainv_min*h*num_ins; % Ainv_min*k*h*Iidf(lmax_cur);

        if( Arep_rem - Arep_ins > 0)
            % If we have space to insert all the repeaters for these XCs, then
            % go ahead and do it
            
            new_xc_length = round(lmax_cur/k);
            Iidf(new_xc_length) = Iidf(new_xc_length) + k*Iidf(lmax_cur);
            Iidf(lmax_cur) = 0;
            
            k_vec(lmax_cur) = k;
            h_vec(lmax_cur) = h;
            lmax_cur = lmax_cur - 1;
            Arep_rem = Arep_rem - Arep_ins;
            Arep_used = Arep_used + Arep_ins;
            num_vec(xc_tier) = num_vec(xc_tier) + num_ins;
            size_vec(xc_tier) = h;
        else
            add_repeaters = 0;
        end
    else
        add_repeaters = 0;
    end
end

%repstr = sprintf('k %d\t h %d\t Arep_ins %d\t Arep_rem %d\t lmax_cur %d \t Iidf(l) %d',k,h,Arep_ins,Arep_rem,lmax_cur,Iidf(lmax_cur));
%disp(repstr)
%% pack outputs
% wire.via_area = A_vias_wiring + A_vias_repeaters;
% wire.via_area_wires = A_vias_wiring;
% wire.via_area_repeaters = A_vias_repeaters;
% wire.area_per_layer = A_layer;
% wire.delay_rc = tau_rc_vec;
% wire.delay_repeaters = tau_rep_vec;

repeater.num_per_wire = k_vec;
repeater.size = h_vec;
repeater.area_total = Arep_used;
repeater.num_per_tier = num_vec;
        
    
    
    