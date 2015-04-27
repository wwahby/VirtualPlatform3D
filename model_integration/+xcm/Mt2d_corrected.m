function [Mt2dc, term3, term4, term4_alt, h, g] = Mt2d_corrected(Lx, Nuc_1d, w_tsv)


Mt2dj = xcm.Mt_2d_joyner(Lx);

h = xcm.calc_h(Lx, Nuc_1d, w_tsv);
g = xcm.calc_g(Lx, Nuc_1d, w_tsv);

% For completely symmetric unit cells, g and h are identical
%term2 = conv(g,g);
term3 = conv(h,h);

N_tsvs = Nuc_1d^2;
term4 = zeros(1,length(Mt2dj));
term4(1:2*w_tsv+1) = N_tsvs*xcm.Mt_2d_joyner(w_tsv);

dfxc_intra = N_tsvs*xcm.Mt_2d_joyner(w_tsv);
dfxc_extra = xcm.Mt_2d_joyner(w_tsv*Nuc_1d);
dfxc_extra(1:2*w_tsv+1) = (1:2*w_tsv+1) - dfxc_intra; % removed within-FZ connections

unit_cell_pitch = Lx/Nuc_1d;
term4_alt = zeros(1,length(Mt2dj));
term4_alt(1:2*w_tsv+1) = dfxc_intra;
for dfxc_ind = 1:length(dfxc_extra)
    dfxc_l = dfxc_ind -1;
    dfxc_m = floor(dfxc_l/unit_cell_pitch) + 1;
    act_ind = dfxc_l + dfxc_m*unit_cell_pitch;
    term4_alt(act_ind) = term4_alt(act_ind) + dfxc_extra(dfxc_ind);
end

Mt2dc = Mt2dj - 2*term3 + term4;