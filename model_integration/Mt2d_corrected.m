function Mt2dc = Mt2d_corrected(Lx, Nuc_1d, w_tsv)


Mt2dj = Mt_2d_joyner(Lx);

h = calc_h(Lx, Nuc_1d, w_tsv);

% For completely symmetric unit cells, g and h are identical
%term2 = conv(g,g);
term3 = conv(h,h);

N_tsvs = Nuc_1d^2;
term4 = zeros(1,length(Mt2dj));
term4(1:2*w_tsv+1) = N_tsvs*Mt_2d_joyner(w_tsv);

Mt2dc = Mt2dj - 2*term3 + term4;