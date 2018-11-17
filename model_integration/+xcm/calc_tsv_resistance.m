function Rtsv = calc_tsv_resistance(rho_m, d_tsv, h)
	Rtsv = rho_m*h/pi/(d_tsv/2)^2;
