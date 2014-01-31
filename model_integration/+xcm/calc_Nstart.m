function Nstart = calc_Nstart(Lx,S,r,g_tsv)

[lmax_3d l3d Ns Ng] = xcm.get_params_3d(Lx,S,r);


Nnst = xcm.calc_Nnst(Lx,S,r,g_tsv);
Nnsb = xcm.calc_Nnsb(Lx,S,r,g_tsv);

Nstart = Ng - Nnst - Nnsb;