%clear all
%close all

% Ng = 16e6;
% S = 4;
% r = 1;

% Ng = 200^2;
% Nuc_1d = 10;
% w_tsv = 5;

Ng = 100^2;
Nuc_1d = 10;
w_tsv = 6;

S = 1;
r = 1;
Ns = Ng/S;
Lx = round(sqrt(Ns));



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
[Mt2dc, term3, term4, term4_alt, h, g, term3bf_h] = xcm.Mt2d_corrected(Lx, Nuc_1d, w_tsv);
Mt2dc_alt = Mt2d - 2*term3 + term4_alt;
Nstart = xcm.calc_Nstart(Lx,S,r,g_tsv);

%h = xcm.calc_h(Lx, Nuc_1d, w_tsv);

%% Brute Force
Mt2dbf = xcm.Mt_2d_brute_force(Lx);
[Mt2dbfc, sfxc_bfc, dfxc_bfc] = xcm.Mt2d_brute_force_corrected(Lx, Nuc_1d, w_tsv);

%% More brute force
[term2bf, term3bf] = xcm.calc_mid_terms_brute_force(Lx, Nuc_1d, w_tsv);
[Mt2dc_bf_mult, Mt2d_bf_mult, term2bf_mult, term3bf_mult] = xcm.calc_Mt2d_brute_force_mult(Lx, Nuc_1d, w_tsv);
%% Plots
Mt2dc_constructed = Mt2d -2*term2bf + term4_alt;

htest = [h zeros(1,length(Mt2d)-length(h))];
term2_alt = 2*term3 - 2*w_tsv*Nuc_1d*htest;
term2_alt(1) = term2_alt(1) + w_tsv^2*Nuc_1d^2;
cor_con = 2*term2_alt - term4_alt;
cor_con = 4*term2bf;
Mt2dc_constructed = Mt2d - cor_con;
% 
% figure(1)
% clf
% loglog(iidf);
% hold on
% %loglog(Iidf_joyner,'r--')
% title('Iidf')
% 
% figure(2)
% clf
% plot(Mt3d,'b')
% hold on
% plot(Mt2d,'k:')
% %plot(Mt_joyner,'r--')
% 
% figure(2)
% clf
% plot(nnst,'b')
% hold on
% %plot(Nnst_joyner,'r--')
% set(gca,'yscale','log')
% title('Nnst')
% 
% figure(3)
% clf
% plot(nnsb,'b')
% hold on
% %plot(Nnsb_joyner,'r--')
% %set(gca,'yscale','log')
% title('Nnsb')
% 
% figure(4)
% clf
% plot(Nstart,'b')
% % plot(Nstart_joyner,'r--')

%%
% figure(5)
% clf
% hold on
% plot(Mt2d,'b-')
% plot(Mt2dc,'r-')
% plot(Mt2dc_alt,'m-')
% plot(Mt2dbf,'g-.')
% plot(Mt2dbfc,'c-')
% plot(Mt2dc_constructed,'r--')
% grid on
% legend('2D','2DC','2DC ALT','2D BF', '2DC BF', '2DC BF2 4ALT','location','s')
% xlabel('Separation (GP)')
% ylabel('Site Function')
% set(gca,'yscale','log')
% set(gca,'xscale','log')
% fixfigs(5,2,14,12)
% 
% figure(6)
% clf
% hold on
% plot(abs(Mt2dbf - Mt2d),'k')
% plot(abs(Mt2dbf - Mt2dc),'r')
% plot(abs(Mt2dbfc - Mt2dc_alt),'m')
% plot(abs(Mt2dbfc - Mt2dc),'g--')
% plot(abs(Mt2dbf - Mt2dbfc),'b')
% grid on
% legend('2DBF vs 2D','2DBF vs 2DC','2DBFC vs 2DC\_ALT','2DBFC vs 2DC','2DBF vs 2DBFC','location','w')
% xlabel('Separation (GP)')
% ylabel('Raw Error in Site Function')
% set(gca,'yscale','log')
% set(gca,'xscale','log')
% fixfigs(6,2,14,12)
% 
figure(7)
clf
hold on
plot(abs(Mt2dbf - Mt2d)./Mt2dbf,'k')
plot(abs(Mt2dbf - Mt2dc)./Mt2dbf,'r')
plot(abs(Mt2dbfc - Mt2dc_alt)./Mt2dbfc,'m')
plot(abs(Mt2dbfc - Mt2dc)./Mt2dbfc,'g-')
plot(abs(Mt2dbf - Mt2dbfc)./Mt2dbf,'b')
plot(abs(Mt2dbfc - Mt2dc_constructed)./Mt2dbfc,'r--')
grid on
legend('2DBF vs 2D','2DBF vs 2DC','2DBFC vs 2DC\_ALT','2DBFC vs 2DC','2DBF vs 2DBFC','location','n')
xlabel('Separation (GP)')
ylabel('Relative Deviation in Site Function')
set(gca,'yscale','log')
set(gca,'xscale','log')
fixfigs(7,2,14,12)
% 
% figure(8)
% clf
% hold on
% plot(abs(2*term3-sfxc_bfc)./sfxc_bfc,'b')
% plot(abs(term4_alt-dfxc_bfc)./dfxc_bfc,'g')
% plot(abs(term4-dfxc_bfc)./dfxc_bfc,'r--')
% xlabel('Separation (GP)')
% ylabel('Relative Deviation in Site Function Corrections')
% set(gca,'yscale','log')
% set(gca,'xscale','log')
% legend('SFXC','DFXC\_ALT','DFXC')
% fixfigs(8,2,14,12)
% 
% figure(9)
% clf
% hold on
% plot(sfxc_bfc,'b')
% plot(2*term3,'b--')
% plot(dfxc_bfc,'r')
% plot(term4_alt,'k')
% plot(term4,'g--')
% plot(sfxc_bfc + dfxc_bfc,'m')
% plot(2*term2bf,'b-.')
% xlabel('Separation (GP)')
% ylabel('Site Function Corrections')
% set(gca,'yscale','log')
% set(gca,'xscale','log')
% legend('SFXC','Term3','DFXC','Term4\_Alt','Term4','SFXC+DFXC','T2+3 BF','location','s')
% grid on
% fixfigs(9,2,14,12)

%%
figure(10)
clf
hold on
plot(term3,'b')
plot(term3bf_h,'r--')
plot(term2bf,'g:')
plot(term3bf,'m--')
 set(gca,'yscale','log')
% set(gca,'xscale','log')
xlabel('Separation (GP)')
ylabel('Site Function Correction')
legend('term3 - conv','term3 - BFH','term2 - BF','term3 - BF')
fixfigs(10,2,14,12)

figure(11)
clf
hold on
plot(Mt2dbf,'color',[0.2 0.2 0.2])
plot(Mt2dbfc,'k')
plot(Mt2d,'b--')
plot(Mt2dc,'b-')
plot(Mt2dc_alt,'r')
plot(Mt2dc_constructed,'g--')
plot(Mt2dc_bf_mult,'color',[1 0.5 0])
plot(Mt2d_bf_mult,'color',[1 0.5 0],'linestyle','--')
xlabel('XC Length (GP)')
ylabel('Site function')
set(gca,'yscale','log')
set(gca,'xscale','log')
legend('2DBF','2DBFC','2D','2DC','2DCA','2DCC','2DCMLT','2DMLT','location','s')
grid on
fixfigs(11,2,14,12)

t3_4x = term3;
t3_4x(2:end) = 2*t3_4x(2:end);
figure(12)
clf
hold on
plot(sfxc_bfc,'b')
plot(dfxc_bfc,'r')
plot(sfxc_bfc + dfxc_bfc,'k')
plot(2*t3_4x,'m')
plot(2*t3_4x - term4_alt,'color',[1 0.6 0])
xlabel('Separation (GP)')
ylabel('Site Function Corrections')
set(gca,'yscale','log')
set(gca,'xscale','log')
legend('SFXC','DFXC','SFXC+DFXC','2T3','2T3-T4','location','s')
grid on
fixfigs(12,2,14,12)

figure(13)
clf
hold on
plot(abs((2*t3_4x-term4_alt)-(sfxc_bfc+dfxc_bfc))./(sfxc_bfc+dfxc_bfc),'b')
plot(abs((2*t3_4x)-(sfxc_bfc+dfxc_bfc))./(sfxc_bfc+dfxc_bfc),'r')
plot(abs((cor_con)-(sfxc_bfc+dfxc_bfc))./(sfxc_bfc+dfxc_bfc),'g')
xlabel('Separation (GP)')
ylabel('Relative Deviation in Site Function Corrections')
%set(gca,'yscale','log')
%set(gca,'xscale','log')
%legend('SFXC','DFXC\_ALT','DFXC')
fixfigs(13,2,14,12)
