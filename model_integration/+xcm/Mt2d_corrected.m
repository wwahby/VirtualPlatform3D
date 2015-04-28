function [Mt2dc, term3, term4, term4_alt, h, g, term3bf] = Mt2d_corrected(Lx, Nuc_1d, w_tsv)


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
dfxc_extra(1:2*w_tsv+1) = dfxc_extra(1:2*w_tsv+1) - dfxc_intra; % removed within-FZ connections

unit_cell_pitch = Lx/Nuc_1d;
term4_alt = zeros(1,length(Mt2dj));
term4_alt(1:2*w_tsv+1) = dfxc_intra;
for dfxc_ind = 1:length(dfxc_extra)
    dfxc_l = dfxc_ind -1;
    dfxc_m = floor(dfxc_l/w_tsv) + 0;
    act_l = dfxc_l + (dfxc_m-1)*unit_cell_pitch;
    act_ind = act_l + 1;
    
    dfxc_l_res = dfxc_l - dfxc_m*w_tsv;
    act_l = dfxc_l_res + dfxc_m * unit_cell_pitch;
    act_ind = act_l + 1;
    term4_alt(act_ind) = term4_alt(act_ind) + dfxc_extra(dfxc_ind);
end

lmax = 2*Lx;
T = Lx/Nuc_1d;
t = 1/2*(T-w_tsv);

term3bf = zeros(1,2*Lx+1);
for l_ind = 1:lmax+1
    l = l_ind -1;
    for lx_ind = 1:l_ind
        lx = lx_ind - 1;
        hbfx = xcm.calc_h_brute_force(lx,t,w_tsv,T,Nuc_1d,Lx);
        hbfy = xcm.calc_h_brute_force(l-lx,t,w_tsv,T,Nuc_1d,Lx);
        term3bf(l_ind) = term3bf(l_ind) + hbfx*hbfy;
    end
end
       
Mt2dc = Mt2dj - 2*term3 + term4;