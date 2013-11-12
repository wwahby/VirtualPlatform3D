function Nc = calc_Nc(Mt,Lx,S,r,g_tsv)

Nstart = calc_Nstart(Lx,S,r,g_tsv);
Nc = Mt./Nstart;