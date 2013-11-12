function Iexp = calc_Iexp(alpha,k,p,Mt,Lx,S,r,g_tsv)


Nc = calc_Nc(Mt,Lx,S,r,g_tsv);
Nb = calc_Nb(Nc);

Iexp = alpha*k./Nc.*( (1+Nb).^p + (Nb+Nc).^p - (1+Nb+Nc).^p - Nb.^p);