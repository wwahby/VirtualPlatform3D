%Power supply noise functioin
%Li Zheng, 9/21/2012

function psn = psn_fast(Nstrata,layer,RTSV,LTSV,RPKG,LPKG,Cd,Jchip,temperature,acell,R,rpad)

% %%% Testing code%%%%%
% clear
% clc
% close all
% Nstrata = 2; %Number of Strata
% layer = 2; %Point we want to calculate
% RTSV = 0.01;  %Resistance of a TSV
% LTSV = 2.5e-11; %Inductance of a TSV
% RPKG = 0.006; %Resistance of package
% LPKG = 0.5e-9; %Inductance of package
% Cd = 0.0053;%Percentage of chip area for decoupling capaciance, 0.x
% Jchip = 1e6;%Current density, A/cm^2
% temperature = 20; %Chip temperature
% acell = 212e-6; %unit cell side length
% R = 0.17; %segment resistance
% rpad = 4.04e-6; %Equivalent pad radius
% %%%%%%%%%%%%%%%%%%%%%

T0 = 20;                                        % Original temperature
T = temperature;                                % Current temperature
alpha = 3.9e-3;                                 % Temperature coefficient of resistance (for copper)
Js = Jchip*ones(1,Nstrata);                     % Current density, A/cm^2
% ef = decap*ones(1,Nstrata);                     % Percentage of decap
% Oth=0.65*10^(-9)*ones(1,Nstrata);               % Gate oxide thickness
% Cd=3.9*8.854*10^(-12)./Oth.*ef;                 % Capacitance density
Cd = Cd.*ones(1,Nstrata);
acell = acell.*ones(1,Nstrata);                 % Capacitance density
R = R*(1+alpha*(T-T0)).*ones(1,Nstrata);        % Segment Resistance
rpad = rpad.*ones(1,Nstrata);                   % Equivalent pad radius
Lv=zeros(1,Nstrata);                            % Package inductance and via inductance per IO
Rv=zeros(1,Nstrata);                            % Package resistance and via resistance per IO

if Nstrata==1
    Lv(1)= LPKG;
    Rv(1)= RPKG*(1+alpha*(T-T0));
end
if Nstrata>1
    Lv(1)=LPKG;
    Rv(1)=RPKG*(1+alpha*(T-T0));
    for j=2:Nstrata
        Lv(j)= LTSV; %u*lv(j)/2/pi*log(1+2.84*lv(j)/pi/rv(j))+4*0.199*u*lv(j)*log(1+1.0438*lv(j)/sqrt(2*acell(j)^2))-4*0.199*u*lv(j)*log(1+1.0438*lv(j)/Padpitch(j));
        Rv(j)= RTSV*(1+alpha*(T-T0));%0.006;
    end
end
xp=acell(layer);    %point within the unit cell
yp=0;

% Define the start frequency and end frequency

f1=5*10^4;
f2=1*10^10;
Nf=1200;
df=log10(f2/f1)/Nf;

AA=zeros(Nstrata,Nstrata);
BB=zeros(Nstrata,1);

Vp=zeros(1,Nf);
Up=zeros(1,Nf);
AbsVp=zeros(1,Nf);
Deg=zeros(1,Nf);
Db=zeros(1,Nf);
%AbsVspv=zeros(Nf,1);


for j=1:1:Nf
    f(j)=f1*10^(df*(j-1));
    Omg=2*pi*f(j);
    s=Omg*i;
    
    lamda=zeros(1,Nstrata);
    for k=1:Nstrata
        lamda(k)=-2*i*Omg*R(k)*Cd(k);
    end
    
    
    if Nstrata==1
        m=1;
        GG = power_noise.Gr1_vec(rpad(m),0,0,0,lamda(m),acell(m));
        AA(m,1)=1+R(m)/(4*s*Lv(m)+4*Rv(m))*GG;
        BB(m,1)=R(m)/(4*s*Lv(m)+4*Rv(m))*GG*Js(m)/2/s/Cd(m);
    elseif Nstrata>1
        for m=1:Nstrata
            GG = power_noise.Gr1_vec(rpad(m),0,0,0,lamda(m),acell(m));
            if m==1
                AA(m,1)=1+R(m)/(4*s*Lv(m)+4*Rv(m))*GG+R(m)/(4*s*Lv(m+1)+4*Rv(m+1))*GG;
                AA(m,2)=-R(m)/(4*s*Lv(m+1)+4*Rv(m+1))*GG;
                BB(m,1)=R(m)/(4*s*Lv(m)+4*Rv(m))*GG*Js(m)/2/s/Cd(m)+R(m)/(4*s*Lv(m+1)+4*Rv(m+1))*GG*Js(m)/2/s/Cd(m)-R(m)/(4*s*Lv(m+1)+4*Rv(m+1))*GG*Js(m+1)/2/s/Cd(m+1);
            end
            
            if m==Nstrata
                AA(m,Nstrata-1)=-R(m)/(4*s*Lv(m)+4*Rv(m))*GG;
                AA(m,Nstrata)=1+R(m)/(4*s*Lv(m)+4*Rv(m))*GG;
                BB(m,1)=-R(m)/(4*s*Lv(m)+4*Rv(m))*GG*Js(m-1)/2/s/Cd(m-1)+R(m)/(4*s*Lv(m)+4*Rv(m))*GG*Js(m)/2/s/Cd(m);
            end
            
            if m>1 && m<Nstrata
                AA(m,m-1)=-R(m)/(4*s*Lv(m)+4*Rv(m))*GG;
                AA(m,m)=1+R(m)/(4*s*Lv(m)+4*Rv(m))*GG+R(m)/(4*s*Lv(m+1)+4*Rv(m+1))*GG;
                AA(m,m+1)=-R(m)/(4*s*Lv(m+1)+4*Rv(m+1))*GG;
                BB(m,1)=-R(m)/(4*s*Lv(m)+4*Rv(m))*GG*Js(m-1)/2/s/Cd(m-1)+R(m)/(4*s*Lv(m)+4*Rv(m))*GG*Js(m)/2/s/Cd(m)+R(m)/(4*s*Lv(m+1)+4*Rv(m+1))*GG*Js(m)/2/s/Cd(m)-R(m)/(4*s*Lv(m+1)+4*Rv(m+1))*GG*Js(m+1)/2/s/Cd(m+1);
            end
        end
    end
    Upad=AA^(-1)*BB;
    
    GG = power_noise.Gr1_vec(xp,yp,0,0,lamda(layer),acell(layer));
    if Nstrata==1&&layer==1
        Up(j)=R(layer)*GG*(-Upad(layer,1)+Js(layer)/2/s/Cd(layer))/(4*s*Lv(layer)+4*Rv(layer));
    end
    
    if Nstrata>1
        if layer==1
            Up(j)=R(layer)*GG*((-Upad(layer,1)+Js(layer)/2/s/Cd(layer))/(4*s*Lv(layer)+4*Rv(layer))-(Upad(layer,1)-Upad(layer+1,1)-Js(layer)/2/s/Cd(layer)+Js(layer+1)/2/s/Cd(layer+1))/(4*s*Lv(layer+1)+4*Rv(layer+1)));
        end
        
        if layer==Nstrata
            Up(j)=R(layer)*GG*(Upad(layer-1,1)-Upad(layer,1)-Js(layer-1)/2/s/Cd(layer-1)+Js(layer)/2/s/Cd(layer))/(4*s*Lv(layer)+4*Rv(layer));
        end
        
        if layer>1 && layer<Nstrata
            Up(j)=R(layer)*GG*((Upad(layer-1,1)-Upad(layer,1)-Js(layer-1)/2/s/Cd(layer-1)+Js(layer)/2/s/Cd(layer))/(4*s*Lv(layer)+4*Rv(layer))-(Upad(layer,1)-Upad(layer+1,1)-Js(layer)/2/s/Cd(layer)+Js(layer+1)/2/s/Cd(layer+1))/(4*s*Lv(layer+1)+4*Rv(layer+1)));
        end
    end
    
    Vp(j)=Up(j)-Js(layer)/2/s/Cd(layer);
    AbsVp(j)=abs(Vp(j));
    Deg(j)=phase(Vp(j))/pi*180;
    Db(j)=20*log10(AbsVp(j));
end

Db0=Db(1);
Deg0=-180;
f0=0;

%%%%% Inverse FFT%%%%%%
%Bumped Npts up from 1024 to 1024*10 due to resolution issues noticed in 3D Sandy Bridge test case when sweeping decap from 0-0.2. WWAHBY 2014-10-06
% Default is 2^10, seems like Npts should be a power of 2
Npts=1024*10; %[FIX] May need a dynamic way to determine FFT resolution. WWAHBY 2014.10.06
tstop=1000e-9;
tstep=tstop/Npts;
%PWL piecewise linear function
v=zeros(1,Npts);
for k=1:1:Npts
    if k<=100e-9/tstep
        v(k)=0;
    end
    if k>=100e-9/tstep && k<=100.6e-9/tstep
        v(k) = (1/(0.6e-9))*(k-(100e-9/tstep))*tstep;
    end
    if k>=100.6e-9/tstep && k<= 900e-9/tstep
        v(k)=1;
    end
    if k>=900e-9/tstep && k<=900.6e-9/tstep
        v(k) = 1-(1/(0.6e-9))*(k-(900e-9/tstep))*tstep;
    end
    if k>900.6e-9/tstep
        v(k)=0;
    end
end
NFFT=2^nextpow2(length(v));
Y = fft(v,NFFT);%/NFFT;
Fs=1/(tstep);
f1 = Fs/2*linspace(0,1,NFFT/2+1);
fft_sys_phase=interp1(f,Deg,f1,'cubic',90);
fft_sys_mag=interp1(f,10.^(Db/20),f1,'cubic',1e-5);
fft_sys=(fft_sys_mag).*exp(1i*deg2rad(fft_sys_phase));

pval_sys=sum(fft_sys_mag.^2);
output=Y(1:NFFT/2+1).*fft_sys;
%opval=ifft(output,NFFT);
opval=ifft(output,NFFT,'symmetric');
%figure(4)
%plot(0:tstep:(length(opval)-1)*tstep,opval)
max_PSN=abs(min(opval))+ mean(sqrt((opval(1:round(90e-9/tstep))).^2));
psn=max_PSN;
end

% file_name=strcat('3.sp');
% fid=fopen(file_name,'W');
% fprintf(fid,'*Laplace voltage gain \n');
% fprintf(fid,'.tran 0.01n 90n \n');
% fprintf(fid,'.print V(output) \n');
% fprintf(fid,'.option post\n');
% fprintf(fid,'Vs input 0 PWL(0 0 10n 0 10.6n 1 100n 1 100.6 0)\n');
% %fprintf(fid,'Vs input 0 AC 1 sin(0 1 1000000000\n');
% %fprintf(fid,'Rs output 0 1\n');
% 
% fprintf(fid,'\n');
% 
% 
% %fprintf(fid,'.op\n')
% 
% %fprintf(fid,'.print V(n_19_22) \n')
% fprintf(fid,'Esource output 0 Freq input 0 \n');
% fprintf(fid,'+ %d     %d     %d\n',f0,Db0,Deg0);
% % fprintf(fid,'+ delf=10000    maxf=100000000000\n');
% for j=1:Nf
%     fprintf(fid,'+ %d     %d     %d\n',f(j),Db(j),Deg(j));
% end
% 
% fprintf(fid,'.end')
% fclose(fid);
% 
% semilogx(f,AbsVp);
% 
% semilogx(f,Deg);