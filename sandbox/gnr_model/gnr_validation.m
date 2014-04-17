gnr_kumar_length = [1 5 10 15 20]*1e-6; % (m) GNR length
gnr_kumar_width = 10e-9; % (m) GNR width

gnr_kumar_layers = [2 4 6 8];
gnr_kumar_R = [13.66 20.05 27.53 34.86 42.26;
               13.66 19.74 24.73 28.62 32.52;
               13.66 19.74 24.18 26.83 29.40;
               13.66 19.74 24.10 26.29 28.16]*1e3;
           
gnr_my_kumar_R = [21.12 103.6 207.1 310.5 414.0]*1e3;
           
gnr_nishad_R = [8.11 40.41 80.81 121.2 161.6;
               4.13 20.20 40.41 61.61 80.82;
               2.85 13.50 26.95 40.42 53.88;
               2.24 10.15 20.23 30.32 40.42]*1e3;
           
figure(1)
set(gca,'DefaultAxesColorOrder',[0 0 0; 0 0 1; 0 1 0; 1 0 0]);
set(gca,'DefaultAxesLineStyleOrder','-')
clf
plot(gnr_kumar_length*1e6,gnr_kumar_R/1e3)

set(gcf,'DefaultAxesColorOrder',[0 0 0; 0 0 1; 0 1 0; 1 0 0]);
set(gcf,'DefaultAxesLineStyleOrder','--')
hold on
plot(gnr_kumar_length*1e6,gnr_nishad_R/1e3,'--')

fixfigs(1,3,14,12)