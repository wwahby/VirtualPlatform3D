function [C_scaled C_fixed C_fixed_venk cap_const] = calc_cu_wire_capacitance(epsrd,pitch,wire_length,aspect_ratio,width_fraction,vertical_spacing)
% epsrd: (-) Dielectric relative permittivity
% pitch: (m) distance from the starting edge of one wire to the starting
% edge of the next wire
% length: (m) length of the wire
% width_fraction: (-) horizontal space between wires / wire pitch
% aspect ratio: (-) wire height / wire width
% vertical spacing: (m) Distance from the top of one wire to the bottom of
% the next wire up

eps0 = 8.854e-12; % (F/m) Vacuum permittivity
eps = eps0*epsrd; % (F/m) Dielectric permittivity


% capacitance assuming vertical_spacing = pitch (one full pitch)
C_scaled = 2*eps*wire_length*width_fraction*(width_fraction + aspect_ratio/(1 - width_fraction));


% capacitance assuming vertical spacing is fixed at the value input
cap_const = zeros(1,length(pitch));

for ind = 1:length(pitch)
    width = pitch(ind)*width_fraction;
    horiz_space = pitch(ind) - width;
    thickness = aspect_ratio*width;
    vert_space = vertical_spacing;
    cap_const(ind) = calc_capacitance_constant_full(width,horiz_space,thickness,vert_space);
end

C_fixed = 2*eps*wire_length*width_fraction*(width_fraction*pitch/vertical_spacing + aspect_ratio/(1-width_fraction));

C_fixed_venk = cap_const * wire_length*eps;

