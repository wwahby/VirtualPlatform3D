function Mt2d = Mt_2d_joyner(Lx)
% Calculates number of gate pairs separated by length l for a 2D chip or
% layer
% Inputs:
%   Lx (gate pitches) - Length of chip in gate pitches (chip is assumed to
%   be square!)
% Outputs:
%   Mt2d: (1x[lmax + 1] vector) # of gate pairs separated by length l

[lmax l Ns] = get_params_2d(Lx);

Mt2d = zeros(1,length(l));

A = l < Lx;
B = (Lx <= l) & (l < (2*Lx - 1));

Mt2d(A) = 2*Ns*l(A) - 2*Lx*l(A).^2 + 1/3*l(A).^3;
Mt2d(1) = Ns;
Mt2d(B) = 1/3 * (2*Lx - l(B)).^3;


