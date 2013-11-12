function Nnst = calc_Nnst(Lx,S,r,g_tsv)

[lmax_3d l3d Ns Ng] = get_params_3d(Lx,S,r);

Nnst = zeros(1,length(l3d));

A = l3d <= Lx/2;
B = (Lx/2 < l3d)&(l3d <= Lx);
C = (Lx < l3d)&(l3d <= 3/2*Lx);
D = (3/2*Lx < l3d)&(l3d<=Ns);

Nnst(A) = l3d(A);
Nnst(B) = l3d(B) + (l3d(B) - Lx/2 - 1).*(l3d(B)-Lx/2);
Nnst(C) = l3d(C)*Lx - 3/4*Ns + 1/2*Lx;
Nnst(D) = Ns-(2*Lx-l3d(D)).*(2*Lx-l3d(D)-1);


if(g_tsv > 0)
    mnst = Nnst./Ns;
    Nnst = Nnst + g_tsv*(1-2*mnst);
end