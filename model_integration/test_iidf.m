%clear all
%close all

Ng = 16e6;
S = 4;
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



iidf = calc_Iidf(alpha,k,p,Lx,S,r);
%AA = calc_Nstart(Lx,S,r);
nnst = calc_Nnst(Lx,S,r);
nnsb = calc_Nnsb(Lx,S,r);
Mt3d = Mt_3d_joyner(Lx,S,r);
Mt2d = Mt_2d_joyner(Lx);
Nstart = calc_Nstart(Lx,S,r);

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