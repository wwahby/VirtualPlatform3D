%clear all
%close all

Ng = 1.6e9;
S = 4;
r = 100;

Ns = Ng/S;
Lx = round(sqrt(Ns));

% Recalculate these to make sure everything is a nice integer
Ns = Lx^2;
Ng = Ns*S;


p = 0.6; % rent exponent
fo = 4; % avg fanout
alpha = fo/(fo+1); % input terminal fraction
k = 3/alpha; %rent constant

Nuc_1d = 200;
T = Lx/Nuc_1d;
w_tsv = 32;
N_tsvs = Nuc_1d^2;
g_tsv = (Nuc_1d*w_tsv)^2;

%%
iidf_c = calc_Iidf_corrected(alpha,k,p,Lx,S,r,Nuc_1d,w_tsv);
iidf_j = calc_Iidf(alpha,k,p,Lx,S,r);

%AA = calc_Nstart(Lx,S,r);
Mt3dj = Mt_3d_joyner(Lx,S,r);
Mt3dc = Mt_3d_corrected(Lx,S,r,Nuc_1d,w_tsv);
Mt2dj = Mt_2d_joyner(Lx);
Mt2dc = Mt2d_corrected(Lx, Nuc_1d, w_tsv);
Nstart = calc_Nstart(Lx,S,r,g_tsv);
Nnst = calc_Nnst(Lx,S,r,g_tsv);
Nnsb = calc_Nnsb(Lx,S,r,g_tsv);
Nc = calc_Nc(Mt3dc,Lx,S,r,g_tsv);
h = calc_h(Lx, Nuc_1d, w_tsv);
term4 = zeros(1,length(Mt2dj));
term4(1:2*w_tsv+1) = N_tsvs*Mt_2d_joyner(w_tsv);
iexp = calc_Iexp(alpha,k,p,Mt3dc,Lx,S,r,g_tsv);

%%
figure(1)
clf
loglog(iidf_j);
hold on
loglog(iidf_c,'r--')
title('Iidf')
%%
figure(2)
clf
plot(Mt3dj,'b')
hold on
plot(Mt3dc,'r:')
title('Mt3d')

figure(3)
clf
plot(Mt2dj,'b')
hold on
plot(Mt2dc,'r:')
title('Mt2d')

figure(4)
clf
plot(Nnst,'b')
hold on
plot(Nnsb,'r--')
title('Nonstarting gates')
%%
err_raw = abs(iidf_c-iidf_j);
figure(5)
clf
loglog(err_raw)
title('raw error')

err_norm = err_raw./iidf_j;
figure(6)
clf
loglog(err_norm)
title('normalized error')