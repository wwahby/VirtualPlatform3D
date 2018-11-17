function Ctsv = calc_tsv_capacitance(epsr_ox, d_tsv, t_ox, h)
    Na = 1e14;
    ni = 1e10;
    T = 273;
	eps0 = 8.854e-12; % F/m
	eps = epsr_ox*eps0;
	
	r_tsv = d_tsv/2;
	r_ox = r_tsv + t_ox;

	Ctsv = 2*pi*eps*h/log( r_ox/r_tsv );