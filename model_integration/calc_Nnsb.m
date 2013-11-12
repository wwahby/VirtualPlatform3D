function Nnsb = calc_Nnsb(Lx,S,r,g_tsv)

[lmax_3d l3d Ns Ng] = get_params_3d(Lx,S,r);


A = max(0,l3d-Lx-r);
B = min(S-1, floor( (l3d-Lx)/r));
B = max(0,B); % Need to ensure B >= 0 or else compact form of Ngpyr gives incorrect values

C = max(0,l3d-3/2*Lx-r);
D = min(S-1, floor( (l3d-3/2*Lx)/r) );
D = max(0,D); % Need to ensure D >= 0 or else compact form of Ngpyr gives incorrect values

Nnsb = calc_Ngpyr(A,B,r) - 2*calc_Ngpyr(C,D,r);

if(g_tsv > 0)
    mnsb = Nnsb./( (S-1)*Ns );
    Nnsb = Nnsb + g_tsv*(1-2*mnsb);
end