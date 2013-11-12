%Unit cell function
%Li Zheng, 9/24/2012

function [lcell,Ppitch,PGpitch,rpad,Rseg,Cd] = Ucell(area,Npad,Ngrid,padsize,rho,Tseg,Wseg,decap)

%%%%%Test codes%%%%%
% clear
% clc
% close all
% area = 184e-6;        %chip area
% Npad = 32*32;         %total number of power or ground pads
% rho = 1.68*10^-8;     %copper resistivity
% Ngrid = 21*21;        %grid fineness
% padsize = 1;          % Pad size is in terms of segment number
% Tseg = 1e-6;          % Thickness of grid segment
% Wseg = 2e-6;          % width of grid segment
%%%%%%%%%%%%%%%%%%%%

lchip = sqrt(area);                 %chip side length (square chip)
Ppitch = lchip/sqrt(Npad);          %distance between two power pads or two ground pads
lcell = Ppitch/2;
PGpitch = sqrt(lcell^2*2);          %distance between two adjacent power and groudn pads
asegment = Ppitch/sqrt(Ngrid-1);    %length of a segment

if padsize == 1
    alfa = 0.2;
end
if padsize > 1
    alfa = 0.5902;
end
rpad = asegment*padsize*alfa;

Rseg = rho*asegment/(Tseg*Wseg);

Oth=0.65*10^(-9);                   % Gate oxide thickness
Cd=3.9*8.854*10^(-12)/Oth*decap;    % Capacitance density
end
    
