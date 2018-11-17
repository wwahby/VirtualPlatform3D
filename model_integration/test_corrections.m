Nt = 10e6;
S = 4;

w_tsv = 20;
r = 10;
Luc_1d = 90;
Nuc_1d = 100;

area_ratio = (w_tsv/Luc_1d)^2;

Lx = Luc_1d*Nuc_1d;
Ns = Lx^2;



 Mt2dc = xcm.Mt2d_corrected(Lx, Nuc_1d, w_tsv);
 Mt2d = xcm.Mt_2d_joyner(Lx);
 
 err_rel = abs(Mt2d-Mt2dc)./Mt2d;
 
 figure(1)
 clf
 hold on
 plot(err_rel,'k', 'linewidth',1.5)
% set(gca, 'yscale','log')
 
 
 %%
 
 Mt3dc = xcm.Mt_3d_corrected(Lx,S,r,Nuc_1d,w_tsv);
 Mt3d = xcm.Mt_3d_joyner(Lx,S,r);
 
 
 err_rel = abs(Mt3d-Mt3dc)./Mt3d;
 
 figure(2)
 clf
 hold on
 plot(err_rel,'k', 'linewidth',1.5)
% set(gca, 'yscale','log')


%%

alpha = 0.5;
k = 3;
p = 0.6;

Iidf = xcm.calc_Iidf(alpha,k,p,Lx,S,r);
Iidfc = xcm.calc_Iidf_corrected(alpha,k,p,Lx,S,r,Nuc_1d,w_tsv);

err_rel = abs(Iidf-Iidfc)./Iidf;
 
figure(3)
clf
hold on
grid on
plot(err_rel(Iidf>1),'k', 'linewidth',1.5)
%set(gca, 'yscale','log')
xlabel('Interconnect Length (GP)')
ylabel('Relative Error (-)')


%% figure(4)

figure(4)
clf
hold on
grid on
plot(Iidf(Iidf>0), 'k', 'linewidth', 1.5)
%plot(Iidfc(Iidf>1), 'r', 'linewidth', 1.5)
set(gca,'yscale','log')
set(gca,'xscale','log')
xlabel('Interconnect Length (GP)')
ylabel('Number of Interconnects')
