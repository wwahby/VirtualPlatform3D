function [lmax l Ns] = get_params_2d(Lx)
Ns = Lx^2;
lmax = 2*Lx;
l = 0:lmax;