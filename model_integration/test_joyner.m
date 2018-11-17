close all
clear all

alpha = 0.67;
k = 3;
p = 0.6;
Lx = 1000;
S = 4;
r = 1;

Lx = 1000;
Nt = Lx^2;

S=1;
wld1 = xcm.calc_Iidf(alpha, k, p, round(sqrt(Nt/S)), S, r);
%Iidf = calc_Iidf(alpha,k,p,Lx,S,r)

S=2;
wld2 = xcm.calc_Iidf(alpha, k, p, round(sqrt(Nt/S)), S, r);

S=3;
wld3 = xcm.calc_Iidf(alpha, k, p, round(sqrt(Nt/S)), S, r);

S=4;
wld4 = xcm.calc_Iidf(alpha, k, p, round(sqrt(Nt/S)), S, r);


figure(1)
clf()
hold on
plot(wld1, 'k', 'linewidth', 1.5)
plot(wld2, 'b', 'linewidth', 1.5)
plot(wld3, 'r', 'linewidth', 1.5)
plot(wld4, 'g', 'linewidth', 1.5)
grid on
set(gca,'xscale','log')
set(gca,'yscale','log')


%%
alpha = 0.67;
k = 3;
p = 0.75;
Lx = 1000;
S = 4;
r = 1;

Lx = 100000;
Nt = Lx^2;

S=1;
wld1 = xcm.calc_Iidf(alpha, k, p, round(sqrt(Nt/S)), S, r);
%Iidf = calc_Iidf(alpha,k,p,Lx,S,r)

S=2;
wld2 = xcm.calc_Iidf(alpha, k, p, round(sqrt(Nt/S)), S, r);

S=3;
wld3 = xcm.calc_Iidf(alpha, k, p, round(sqrt(Nt/S)), S, r);

S=4;
wld4 = xcm.calc_Iidf(alpha, k, p, round(sqrt(Nt/S)), S, r);


figure(1)
clf()
hold on
plot(wld1, 'k', 'linewidth', 1.5)
plot(wld2, 'b', 'linewidth', 1.5)
plot(wld3, 'r', 'linewidth', 1.5)
plot(wld4, 'g', 'linewidth', 1.5)
grid on
set(gca,'xscale','log')
set(gca,'yscale','log')
