function Iidf = calc_Iidf(alpha,k,p,Lx,S,r)

Mt = Mt_3d_joyner(Lx,S,r);

g_tsv=0;
Iexp = calc_Iexp(alpha,k,p,Mt,Lx,S,r,g_tsv);

Iidf = Mt.*Iexp;