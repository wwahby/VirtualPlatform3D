function [W gate_pitch] = size_gates_nand2(Wmin,Rnand,Cnand,Cint,fo,logic_depth,Tclk,Beta_g)
% Calculates appropriate size for transistors in logic design in order to
% meet a particular clock frequency
% =====
% INPUTS
% =====
% Wmin - (m) Minimum allowable feature size
% Rnand - (Ohm) avg drive resistance of a MINIMUM SIZE nand2 gate
% Cnand - (F) capacitance of a MINIMUM SIZE nand2 gate
% Cint - (F) capacitance of an average length wire\
% fo - (-) avg fanout
% Ld - (-) Logic depth
% Tclk - (s) Clock period
% Beta_g - (-) Fraction of clock period allowed for gate transitions (~0.75)


chi = 4/(fo+3);

W = 0.*logic_depth*fo*chi*Rnand*Cint/(Beta_g*Tclk-0.7*logic_depth*fo*Rnand*Cnand);

if(W < Wmin)
    W = Wmin;
end

gate_pitch = 16*Wmin*sqrt(W/Wmin);