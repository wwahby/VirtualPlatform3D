%clear all
%close all

% Ng = 16e6;
% S = 4;
% r = 1;

Ng = 200^2;
S = 1;
r = 1;
Ns = Ng/S;
Lx = round(sqrt(Ns));

Nuc_1d = 10;
w_tsv = 5;

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
[Mt2dc, term3, term4, term4_alt, h, g, term3bf] = xcm.Mt2d_corrected(Lx, Nuc_1d, w_tsv);
Mt2dc_alt = Mt2d - conv(g,g) - conv(h,h) + term4_alt;
Nstart = xcm.calc_Nstart(Lx,S,r,g_tsv);

%h = xcm.calc_h(Lx, Nuc_1d, w_tsv);

%% Brute Force
Mt2dbf = xcm.Mt_2d_brute_force(Lx);
[Mt2dbfc, sfxc_bfc, dfxc_bfc] = xcm.Mt2d_brute_force_corrected(Lx, Nuc_1d, w_tsv);

%% Plots
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
plot(Mt2dc_alt,'m--')
plot(Mt2dbf,'g-.')
plot(Mt2dbfc,'c--')
grid on
legend('2D','2DC','2DBF','2DBFC')
xlabel('Separation (GP)')
ylabel('Site Function')
set(gca,'yscale','log')
set(gca,'xscale','log')
fixfigs(5,2,14,12)

figure(6)
clf
hold on
plot(abs(Mt2dbf - Mt2d),'k')
plot(abs(Mt2dbf - Mt2dc),'r')
plot(abs(Mt2dbfc - Mt2dc_alt),'m')
plot(abs(Mt2dbfc - Mt2dc),'g--')
plot(abs(Mt2dbf - Mt2dbfc),'b')
grid on
legend('2DBF vs 2D','2DBF vs 2DC','2DBFC vs 2DC\_ALT','2DBFC vs 2DC','2DBF vs 2DBFC','location','w')
xlabel('Separation (GP)')
ylabel('Raw Error in Site Function')
set(gca,'yscale','log')
set(gca,'xscale','log')
fixfigs(6,2,14,12)

figure(7)
clf
hold on
plot(abs(Mt2dbf - Mt2d)./Mt2dbf,'k')
plot(abs(Mt2dbf - Mt2dc)./Mt2dbf,'r')
plot(abs(Mt2dbfc - Mt2dc_alt)./Mt2dbfc,'m')
plot(abs(Mt2dbfc - Mt2dc)./Mt2dbfc,'g--')
plot(abs(Mt2dbf - Mt2dbfc)./Mt2dbf,'b')
grid on
legend('2DBF vs 2D','2DBF vs 2DC','2DBFC vs 2DC\_ALT','2DBFC vs 2DC','2DBF vs 2DBFC','location','n')
xlabel('Separation (GP)')
ylabel('Relative Deviation in Site Function')
set(gca,'yscale','log')
set(gca,'xscale','log')
fixfigs(7,2,14,12)

figure(8)
clf
hold on
plot(abs(2*term3-sfxc_bfc)./sfxc_bfc,'b')
plot(abs(term4_alt-dfxc_bfc)./dfxc_bfc,'g')
plot(abs(term4-dfxc_bfc)./dfxc_bfc,'r--')
xlabel('Separation (GP)')
ylabel('Relative Deviation in Site Function Corrections')
set(gca,'yscale','log')
set(gca,'xscale','log')
legend('SFXC','DFXC\_ALT','DFXC')
fixfigs(8,2,14,12)

figure(9)
clf
hold on
plot(sfxc_bfc,'b')
plot(2*term3,'b--')
plot(dfxc_bfc,'r')
plot(term4_alt,'k')
plot(term4,'g--')
plot(sfxc_bfc + dfxc_bfc,'m')
xlabel('Separation (GP)')
ylabel('Site Function Corrections')
set(gca,'yscale','log')
set(gca,'xscale','log')
legend('SFXC','Term3','DFXC','Term4\_Alt','Term4','SFXC+DFXC','location','s')
fixfigs(9,2,14,12)

figure(10)
clf
hold on
plot(term3,'b')
plot(term3bf,'r--')
% set(gca,'yscale','log')
% set(gca,'xscale','log')
xlabel('Separation (GP)')
ylabel('Site Function Correction')
legend('term3 - conv','term3 - BF')
fixfigs(10,2,14,12)