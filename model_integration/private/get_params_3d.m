function [lmax_3d l3d Ns Ng] = get_params_3d(Lx,S,r)

Ns = Lx^2;
Ng = Ns*S;
lmax_3d = 2*Lx + (S-1)*r;
l3d = 0:lmax_3d;