%clear all
%close all

Ng = 16e6;
S = 4;
r = 1;

Ng = 100;

Ns = Ng/S;
Lx = round(sqrt(Ns));

Nuc_1d = 2;
w_tsv = 1;

% Recalculate these to make sure everything is a nice integer
Ns = Lx^2;
Ng = Ns*S;


p = 0.6; % rent exponent
fo = 4; % avg fanout
alpha = fo/(fo+1); % input terminal fraction
k = 3/alpha; %rent constant


g_tsv = 0;

iidf = xcm.calc_Iidf(alpha,k,p,Lx,S,r);
%AA = calc_Nstart(Lx,S,r);
nnst = xcm.calc_Nnst(Lx,S,r,g_tsv);
nnsb = xcm.calc_Nnsb(Lx,S,r,g_tsv);
Mt3d = xcm.Mt_3d_joyner(Lx,S,r);
Mt2d = xcm.Mt_2d_joyner(Lx);
Mt2dc = xcm.Mt2d_corrected(Lx, Nuc_1d, w_tsv);
Nstart = xcm.calc_Nstart(Lx,S,r,g_tsv);

h = xcm.calc_h(Lx, Nuc_1d, w_tsv);

figure(1)
clf
loglog(iidf);
hold on
%loglog(Iidf_joyner,'r--')
title('Iidf')

figure(2)
clf
plot(Mt3d,'b')
hold on
plot(Mt2d,'k:')
%plot(Mt_joyner,'r--')

figure(2)
clf
plot(nnst,'b')
hold on
%plot(Nnst_joyner,'r--')
set(gca,'yscale','log')
title('Nnst')

figure(3)
clf
plot(nnsb,'b')
hold on
%plot(Nnsb_joyner,'r--')
%set(gca,'yscale','log')
title('Nnsb')

figure(4)
clf
plot(Nstart,'b')
% plot(Nstart_joyner,'r--')


figure(5)
clf
hold on
plot(Mt2d,'b-')
plot(Mt2dc,'r--')
% 
% figure(6)
% clf
% hold on
% plot(h)