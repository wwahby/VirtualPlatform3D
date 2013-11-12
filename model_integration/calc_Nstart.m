function Nstart = calc_Nstart(Lx,S,r,g_tsv)

[lmax_3d l3d Ns Ng] = get_params_3d(Lx,S,r);


Nnst = calc_Nnst(Lx,S,r,g_tsv);
Nnsb = calc_Nnsb(Lx,S,r,g_tsv);

Nstart = Ng - Nnst - Nnsb;