function g = calc_g(Lx, Nuc_1d, w_tsv)
% Lx (GP) Length of chip in gate pitches
% Nuc_1d (-) Number of unit cells (each containing one TSV) along the
%              length of the chip
% w_tsv (GP) Width of the tsvs in gate pitches

% Don't bother doing this if we don't have any TSVs
if(Nuc_1d > 0)
    l = 0:Lx;
    T = Lx/Nuc_1d;
    t = 1/2*(T-w_tsv);

    g = zeros(1,length(l));
    lxm = mod(l,T);
    nt = Nuc_1d - floor(l/T);

    A = lxm < T-t-w_tsv;
    B = lxm > T-t;
    C = ~(A | B);

    g(A) = nt(A)*w_tsv;
    g(B) = (nt(B)-1)*w_tsv;
    g(C) = nt(C)*w_tsv - (lxm(C) -(T-t-w_tsv));
else
    g = zeros(1,Lx+1);
end


% function h = calc_h(lx,t,w,T,Ntsvs,lxmax)
% % lx - 1d tsv length
% % t - TSV offset within cell (in gate lengths)
% % w - TSV width (in gate lengths)
% % T - cell period (in gate lengths)
% % Ntsvs - number of tsvs in 1d
% 
% lmt = mod(lx,T);
% nt = Ntsvs-floor(lx/T);
% 
% if (lx < 0) || (lx > lxmax)
%     h = 0;
% elseif lmt < t
%     h = nt*w;
% elseif lmt > t+w
%     h = (nt-1)*w;
% else
%     h = nt*w-(lmt-t);
% end