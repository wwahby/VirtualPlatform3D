function Mt3d = Mt_3d_corrected(Lx,S,r,Nuc_1d,w_tsv)
% Calculates number of gates separated by length l, in 3D
% Inputs
%   Lx: (GP) length of one side of the chip
%   S: (-) number of device layers
%   r: (GP) Interlayer separation, in gate pitches
%
% Outputs
%   Mt3d (1xlmax+1 vector)
[lmax_2d l2d Ns] = get_params_2d(Lx);

[lmax3d l3d Ns Ng] = get_params_3d(Lx,S,r);

% Get 2D distribution
%Mt2d = Mt_2d_joyner(Lx);
Mt2d = Mt2d_corrected(Lx, Nuc_1d, w_tsv);

Mt3d = zeros(1,length(l3d));
% v=0 case
Mt3d(1:length(l2d)) = (1-dd(l2d)).*S.*Mt2d;

% The rest of the cases
Sind_max = S-1;
for v=1:Sind_max
    MM = [zeros(1,v*r) Mt2d zeros(1,(Sind_max-v)*r) ];
    Mt3d= Mt3d + (2-dd(l3d-v*r)) .* (S-v) .* MM;
end