function Iidf = calc_Iidf_corrected(alpha,k,p,Lx,S,r,Nuc_1d,w_tsv)

Mt = xcm.Mt_3d_corrected(Lx,S,r,Nuc_1d,w_tsv);

g_tsv = (Nuc_1d*w_tsv)^2;
Iexp = xcm.calc_Iexp(alpha,k,p,Mt,Lx,S,r,g_tsv);

Iidf = Mt.*Iexp;
Iidf(isnan(Iidf)) = 0;