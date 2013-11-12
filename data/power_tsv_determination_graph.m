npads = [1024 , 1764 , 3025 , 5184 , 8836 , 7225 , 5929 , 4761 , 5184 , 5776 , 6400 ];
psnm = [3.15E-01 , 1.83E-01 , 1.11E-01 , 7.96E-02 , 5.34E-02 , 6.44E-02 , 7.39E-02 , 8.30E-02 , 7.96E-02 , 7.51E-02 , 7.03E-02 ];
target = 0.075*ones(1,length(psnm));

figure(1)
clf
figure(1)
plot(npads,'b')
xlabel('iteration')
ylabel('Power TSVs Required')

figure(2)
clf
plot(psnm*1e3,'b')
hold on
plot(target*1e3,'k--')
xlabel('iteration')
ylabel('Power supply noise (mV)')
axis tight
ylim([0 350])
xlim([0 11])
xlim([1 11])

fixfigs(1:2,3,14,12)