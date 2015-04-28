function h = calc_h_brute_force(lx,t,w,T,Ntsvs,lxmax)
% lx - 1d tsv length
% t - TSV offset within cell (in gate lengths)
% w - TSV width (in gate lengths)
% T - cell period (in gate lengths)
% Ntsvs - number of tsvs in 1d

lmt = mod(lx,T);
nt = Ntsvs-floor(lx/T);

if (lx < 0) || (lx > lxmax)
    h = 0;
elseif lmt < t
    h = nt*w;
elseif lmt > t+w
    h = (nt-1)*w;
else
    h = nt*w-(lmt-t);
end