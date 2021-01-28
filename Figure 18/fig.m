clearvars
clc
clear -all

%% dt = -212
load('FT_LFS3Hz_1200p_SKF100m_m212.mat')
load('FT_LFS3Hz_1200p_SKF100m_m212_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF100m_m212.mat')
load('fEPSP_LFS3Hz_1200p_SKF100m_m212_keinDA.mat')

fEPSP_1 = fEPSP_LFS3Hz_1200p_SKF100m_m212;
Ft_1 = FT_LFS3Hz_1200p_SKF100m_m212;
fEPSP_keinDA_1 = fEPSP_LFS3Hz_1200p_SKF100m_m212_keinDA;
Ft_keinDA_1 = FT_LFS3Hz_1200p_SKF100m_m212_keinDA;          


load('FT_LFS3Hz_1200p_SKF50m_m212.mat')
load('FT_LFS3Hz_1200p_SKF50m_m212_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF50m_m212.mat')
load('fEPSP_LFS3Hz_1200p_SKF50m_m212_keinDA.mat')

fEPSP_2 = fEPSP_LFS3Hz_1200p_SKF50m_m212;
Ft_2 = FT_LFS3Hz_1200p_SKF50m_m212;
fEPSP_keinDA_2 = fEPSP_LFS3Hz_1200p_SKF50m_m212_keinDA;
Ft_keinDA_2 = FT_LFS3Hz_1200p_SKF50m_m212_keinDA;      


load('FT_LFS3Hz_1200p_SKF25m_m212.mat')
load('FT_LFS3Hz_1200p_SKF25m_m212_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF25m_m212.mat')
load('fEPSP_LFS3Hz_1200p_SKF25m_m212_keinDA.mat')

fEPSP_3 = fEPSP_LFS3Hz_1200p_SKF25m_m212;
Ft_3 = FT_LFS3Hz_1200p_SKF25m_m212;
fEPSP_keinDA_3 = fEPSP_LFS3Hz_1200p_SKF25m_m212_keinDA;
Ft_keinDA_3 = FT_LFS3Hz_1200p_SKF25m_m212_keinDA;      

load('FT_LFS3Hz_1200p_SKF15m_m212.mat')
load('FT_LFS3Hz_1200p_SKF15m_m212_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF15m_m212.mat')
load('fEPSP_LFS3Hz_1200p_SKF15m_m212_keinDA.mat')

fEPSP_4 = fEPSP_LFS3Hz_1200p_SKF15m_m212;
Ft_4 = FT_LFS3Hz_1200p_SKF15m_m212;
fEPSP_keinDA_4 = fEPSP_LFS3Hz_1200p_SKF15m_m212_keinDA;
Ft_keinDA_4 = FT_LFS3Hz_1200p_SKF15m_m212_keinDA;      


load('FT_LFS3Hz_1200p_SKF10m_m212.mat')
load('FT_LFS3Hz_1200p_SKF10m_m212_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF10m_m212.mat')
load('fEPSP_LFS3Hz_1200p_SKF10m_m212_keinDA.mat')

fEPSP_5 = fEPSP_LFS3Hz_1200p_SKF10m_m212;
Ft_5 = FT_LFS3Hz_1200p_SKF10m_m212;
fEPSP_keinDA_5 = fEPSP_LFS3Hz_1200p_SKF10m_m212_keinDA;
Ft_keinDA_5 = FT_LFS3Hz_1200p_SKF10m_m212_keinDA;      


load('FT_LFS3Hz_1200p_SKF5m_m212.mat')
load('FT_LFS3Hz_1200p_SKF5m_m212_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF5m_m212_m0.mat')
load('fEPSP_LFS3Hz_1200p_SKF5m_m212_keinDA.mat')

fEPSP_6 = fEPSP_LFS3Hz_1200p_SKF5m_m212;
Ft_6 = FT_LFS3Hz_1200p_SKF5m_m212;
fEPSP_keinDA_6 = fEPSP_LFS3Hz_1200p_SKF5m_m212_keinDA;
Ft_keinDA_6 = FT_LFS3Hz_1200p_SKF5m_m212_keinDA;      

load('FT_LFS3Hz_1200p_SKF2m_m212.mat')
load('FT_LFS3Hz_1200p_SKF2m_m212_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF2m_m212_m0.mat')
load('fEPSP_LFS3Hz_1200p_SKF2m_m212_keinDA.mat')

fEPSP_7 = fEPSP_LFS3Hz_1200p_SKF2m_m212;
Ft_7 = FT_LFS3Hz_1200p_SKF2m_m212;
fEPSP_keinDA_7 = fEPSP_LFS3Hz_1200p_SKF2m_m212_keinDA;
Ft_keinDA_7 = FT_LFS3Hz_1200p_SKF2m_m212_keinDA;      

load('FT_LFS3Hz_1200p_SKF1m_m212.mat')
load('FT_LFS3Hz_1200p_SKF1m_m212_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF1m_m212_m0.mat')
load('fEPSP_LFS3Hz_1200p_SKF1m_m212_keinDA.mat')

fEPSP_8 = fEPSP_LFS3Hz_1200p_SKF1m_m212;
Ft_8 = FT_LFS3Hz_1200p_SKF1m_m212;
fEPSP_keinDA_8 = fEPSP_LFS3Hz_1200p_SKF1m_m212_keinDA;
Ft_keinDA_8 = FT_LFS3Hz_1200p_SKF1m_m212_keinDA;  

%% dt = -30
load('FT_LFS3Hz_1200p_SKF100m_m30.mat')
load('FT_LFS3Hz_1200p_SKF100m_m30_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF100m_m30.mat')
load('fEPSP_LFS3Hz_1200p_SKF100m_m30_keinDA.mat')

fEPSP_21 = fEPSP_LFS3Hz_1200p_SKF100m_m30;
Ft_21 = FT_LFS3Hz_1200p_SKF100m_m30;
fEPSP_keinDA_21 = fEPSP_LFS3Hz_1200p_SKF100m_m30_keinDA;
Ft_keinDA_21 = FT_LFS3Hz_1200p_SKF100m_m30_keinDA;          


load('FT_LFS3Hz_1200p_SKF50m_m30.mat')
load('FT_LFS3Hz_1200p_SKF50m_m30_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF50m_m30.mat')
load('fEPSP_LFS3Hz_1200p_SKF50m_m30_keinDA.mat')

fEPSP_22 = fEPSP_LFS3Hz_1200p_SKF50m_m30;
Ft_22 = FT_LFS3Hz_1200p_SKF50m_m30;
fEPSP_keinDA_22 = fEPSP_LFS3Hz_1200p_SKF50m_m30_keinDA;
Ft_keinDA_22 = FT_LFS3Hz_1200p_SKF50m_m30_keinDA;      


load('FT_LFS3Hz_1200p_SKF25m_m30.mat')
load('FT_LFS3Hz_1200p_SKF25m_m30_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF25m_m30.mat')
load('fEPSP_LFS3Hz_1200p_SKF25m_m30_keinDA.mat')

fEPSP_23 = fEPSP_LFS3Hz_1200p_SKF25m_m30;
Ft_23 = FT_LFS3Hz_1200p_SKF25m_m30;
fEPSP_keinDA_23 = fEPSP_LFS3Hz_1200p_SKF25m_m30_keinDA;
Ft_keinDA_23 = FT_LFS3Hz_1200p_SKF25m_m30_keinDA;      

load('FT_LFS3Hz_1200p_SKF15m_m30.mat')
load('FT_LFS3Hz_1200p_SKF15m_m230_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF15m_m30.mat')
load('fEPSP_LFS3Hz_1200p_SKF15m_m30_keinDA.mat')

fEPSP_24 = fEPSP_LFS3Hz_1200p_SKF15m_m30;
Ft_24 = FT_LFS3Hz_1200p_SKF15m_m30;
fEPSP_keinDA_24 = fEPSP_LFS3Hz_1200p_SKF15m_m30_keinDA;
Ft_keinDA_24 = FT_LFS3Hz_1200p_SKF15m_m30_keinDA;      


load('FT_LFS3Hz_1200p_SKF10m_m30.mat')
load('FT_LFS3Hz_1200p_SKF10m_m30_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF10m_m30.mat')
load('fEPSP_LFS3Hz_1200p_SKF10m_m30_keinDA.mat')

fEPSP_25 = fEPSP_LFS3Hz_1200p_SKF10m_m30;
Ft_25 = FT_LFS3Hz_1200p_SKF10m_m30;
fEPSP_keinDA_25 = fEPSP_LFS3Hz_1200p_SKF10m_m30_keinDA;
Ft_keinDA_25 = FT_LFS3Hz_1200p_SKF10m_m30_keinDA;      


load('FT_LFS3Hz_1200p_SKF5m_m30.mat')
load('FT_LFS3Hz_1200p_SKF5m_m30_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF5m_m30.mat')
load('fEPSP_LFS3Hz_1200p_SKF5m_m30_keinDA.mat')

fEPSP_26 = fEPSP_LFS3Hz_1200p_SKF5m_m30;
Ft_26 = FT_LFS3Hz_1200p_SKF5m_m30;
fEPSP_keinDA_26 = fEPSP_LFS3Hz_1200p_SKF5m_m30_keinDA;
Ft_keinDA_26 = FT_LFS3Hz_1200p_SKF5m_m30_keinDA;      

load('FT_LFS3Hz_1200p_SKF2m_m30.mat')
load('FT_LFS3Hz_1200p_SKF2m_m30_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF2m_m30.mat')
load('fEPSP_LFS3Hz_1200p_SKF2m_m30_keinDA.mat')

fEPSP_27 = fEPSP_LFS3Hz_1200p_SKF2m_m30;
Ft_27 = FT_LFS3Hz_1200p_SKF2m_m30;
fEPSP_keinDA_27 = fEPSP_LFS3Hz_1200p_SKF2m_m30_keinDA;
Ft_keinDA_27 = FT_LFS3Hz_1200p_SKF2m_m30_keinDA;      

load('FT_LFS3Hz_1200p_SKF1m_m30.mat')
load('FT_LFS3Hz_1200p_SKF1m_m30_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF1m_m30.mat')
load('fEPSP_LFS3Hz_1200p_SKF1m_m30_keinDA.mat')

fEPSP_28 = fEPSP_LFS3Hz_1200p_SKF1m_m30;
Ft_28 = FT_LFS3Hz_1200p_SKF1m_m30;
fEPSP_keinDA_28 = fEPSP_LFS3Hz_1200p_SKF1m_m30_keinDA;
Ft_keinDA_28 = FT_LFS3Hz_1200p_SKF1m_m30_keinDA;  

%% dt = -20
load('FT_LFS3Hz_1200p_SKF100m_m20.mat')
load('FT_LFS3Hz_1200p_SKF100m_m20_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF100m_m20.mat')
load('fEPSP_LFS3Hz_1200p_SKF100m_m20_keinDA.mat')

fEPSP_31 = fEPSP_LFS3Hz_1200p_SKF100m_m20;
Ft_31 = FT_LFS3Hz_1200p_SKF100m_m20;
fEPSP_keinDA_31 = fEPSP_LFS3Hz_1200p_SKF100m_m20_keinDA;
Ft_keinDA_31 = FT_LFS3Hz_1200p_SKF100m_m20_keinDA;          


load('FT_LFS3Hz_1200p_SKF50m_m20.mat')
load('FT_LFS3Hz_1200p_SKF50m_m20_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF50m_m20.mat')
load('fEPSP_LFS3Hz_1200p_SKF50m_m20_keinDA.mat')

fEPSP_32 = fEPSP_LFS3Hz_1200p_SKF50m_m20;
Ft_32 = FT_LFS3Hz_1200p_SKF50m_m20;
fEPSP_keinDA_32 = fEPSP_LFS3Hz_1200p_SKF50m_m20_keinDA;
Ft_keinDA_32 = FT_LFS3Hz_1200p_SKF50m_m20_keinDA;      


load('FT_LFS3Hz_1200p_SKF25m_m20.mat')
load('FT_LFS3Hz_1200p_SKF25m_m20_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF25m_m20.mat')
load('fEPSP_LFS3Hz_1200p_SKF25m_m20_keinDA.mat')

fEPSP_33 = fEPSP_LFS3Hz_1200p_SKF25m_m20;
Ft_33 = FT_LFS3Hz_1200p_SKF25m_m20;
fEPSP_keinDA_33 = fEPSP_LFS3Hz_1200p_SKF25m_m20_keinDA;
Ft_keinDA_33 = FT_LFS3Hz_1200p_SKF25m_m20_keinDA;      

load('FT_LFS3Hz_1200p_SKF15m_m20.mat')
load('FT_LFS3Hz_1200p_SKF15m_m20_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF15m_m20.mat')
load('fEPSP_LFS3Hz_1200p_SKF15m_m20_keinDA.mat')

fEPSP_34 = fEPSP_LFS3Hz_1200p_SKF15m_m20;
Ft_34 = FT_LFS3Hz_1200p_SKF15m_m20;
fEPSP_keinDA_34 = fEPSP_LFS3Hz_1200p_SKF15m_m20_keinDA;
Ft_keinDA_34 = FT_LFS3Hz_1200p_SKF15m_m20_keinDA;      


load('FT_LFS3Hz_1200p_SKF10m_m20.mat')
load('FT_LFS3Hz_1200p_SKF10m_m20_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF10m_m20.mat')
load('fEPSP_LFS3Hz_1200p_SKF10m_m20_keinDA.mat')

fEPSP_35 = fEPSP_LFS3Hz_1200p_SKF10m_m20;
Ft_35 = FT_LFS3Hz_1200p_SKF10m_m20;
fEPSP_keinDA_35 = fEPSP_LFS3Hz_1200p_SKF10m_m20_keinDA;
Ft_keinDA_35 = FT_LFS3Hz_1200p_SKF10m_m20_keinDA;      


load('FT_LFS3Hz_1200p_SKF5m_m20.mat')
load('FT_LFS3Hz_1200p_SKF5m_m20_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF5m_m20.mat')
load('fEPSP_LFS3Hz_1200p_SKF5m_m20_keinDA.mat')

fEPSP_36 = fEPSP_LFS3Hz_1200p_SKF5m_m20;
Ft_36 = FT_LFS3Hz_1200p_SKF5m_m20;
fEPSP_keinDA_36 = fEPSP_LFS3Hz_1200p_SKF5m_m20_keinDA;
Ft_keinDA_36 = FT_LFS3Hz_1200p_SKF5m_m20_keinDA;      

load('FT_LFS3Hz_1200p_SKF2m_m20.mat')
load('FT_LFS3Hz_1200p_SKF2m_m20_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF2m_m20.mat')
load('fEPSP_LFS3Hz_1200p_SKF2m_m20_keinDA.mat')

fEPSP_37 = fEPSP_LFS3Hz_1200p_SKF2m_m20;
Ft_37 = FT_LFS3Hz_1200p_SKF2m_m20;
fEPSP_keinDA_37 = fEPSP_LFS3Hz_1200p_SKF2m_m20_keinDA;
Ft_keinDA_37 = FT_LFS3Hz_1200p_SKF2m_m20_keinDA;      

load('FT_LFS3Hz_1200p_SKF1m_m20.mat')
load('FT_LFS3Hz_1200p_SKF1m_m20_keinDA.mat')
load('fEPSP_LFS3Hz_1200p_SKF1m_m20.mat')
load('fEPSP_LFS3Hz_1200p_SKF1m_m20_keinDA.mat')

fEPSP_38 = fEPSP_LFS3Hz_1200p_SKF1m_m20;
Ft_38 = FT_LFS3Hz_1200p_SKF1m_m20;
fEPSP_keinDA_38 = fEPSP_LFS3Hz_1200p_SKF1m_m20_keinDA;
Ft_keinDA_38 = FT_LFS3Hz_1200p_SKF1m_m20_keinDA;  



ls = 1.5;
ms = 6;

HFS = 215*[1,1]; 
SKF = [3,23];

nopts = 5;

N1 = length(Ft_1);
a = 110;
b = 112;
I1 = [1:nopts:N1];I1 = cat(2,[1:nopts:a],[b:nopts:N1]);


figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
subplot(3,1,1)
plot(Ft_1(I1)/(60*1000),fEPSP_1(I1)-fEPSP_keinDA_1(I1),'ks-',...Ft_2(I1)/(60*1000),fEPSP_2(I1)-fEPSP_keinDA_2(I1),'m^-',...
     Ft_3(I1)/(60*1000),fEPSP_3(I1)-fEPSP_keinDA_3(I1),'ro-',Ft_4(I1)/(60*1000),fEPSP_4(I1)-fEPSP_keinDA_4(I1),'b<-',...
     Ft_5(I1)/(60*1000),fEPSP_5(I1)-fEPSP_keinDA_5(I1),'g*-',Ft_6(I1)/(60*1000),fEPSP_6(I1)-fEPSP_keinDA_6(I1),'c>-',...
     Ft_7(I1)/(60*1000),fEPSP_7(I1)-fEPSP_keinDA_7(I1),'kd-',Ft_8(I1)/(60*1000),fEPSP_8(I1)-fEPSP_keinDA_8(I1),'m^-',...
     HFS,[55,62],'k-',SKF,[25,25],'b',...
     'LineWidth',ls,'MarkerSize',ms)
% plot(Ft_1/(60*1000)-2,fEPSP_1-fEPSP_keinDA_1,'ks-',FT_HFS_test_SKF50m_m200/(60*1000)-1,fEPSP_HFS_test_SKF50m_m200-fEPSP_HFS_test_keinDA_SKF50m_m200,'rs-','LineWidth',ls,'MarkerSize',ms)
legend('100 \muM SKF 38393','25 \muM SKF 38393','15 \muM SKF 38393','10 \muM SKF 38393','5 \muM SKF 38393','2 \muM SKF 38393','1 \muM SKF 38393')
ylabel('\Delta fEPSP')
xlabel('Time (min)')
% title('\Deltat = -212')
xlim([0,410])
ylim([-20,65])
xticks([0:50:450])


%  Ft_1(I1b)/(60*1000),fEPSP_1(I1b)-fEPSP_keinDA_1(I1b),'ks-',Ft_2(I1)/(60*1000),fEPSP_2(I1b)-fEPSP_keinDA_2(I1b),'m^-',Ft_3(I1b)/(60*1000),...
%      fEPSP_3(I1b)-fEPSP_keinDA_3(I1b),'ro-',Ft_4(I1b)/(60*1000),fEPSP_4(I1b)-fEPSP_keinDA_4(I1b),'b<-',...
%      Ft_5(I1b)/(60*1000),fEPSP_5(I1b)-fEPSP_keinDA_5(I1b),'g*-',Ft_6(I1b)/(60*1000),fEPSP_6(I1b)-fEPSP_keinDA_6(I1b),'c>-',...
%      Ft_7(I1b)/(60*1000),fEPSP_7(I1b)-fEPSP_keinDA_7(I1b),'kd-',Ft_8(I1b)/(60*1000),fEPSP_8(I1b)-fEPSP_keinDA_8(I1b),'r*-',...
%%


HFS = 33*[1,1]; 
SKF = [3,23];

N1 = length(Ft_21);
a = 110;
b = 112;
I1 = [1:nopts:N1];
% I1 = cat(2,[1:nopts:a],[b:nopts:N1]);

% figure(2)%'DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
subplot(3,1,2)
plot(Ft_21(I1)/(60*1000),fEPSP_21(I1)-fEPSP_keinDA_21(I1),'ks-',...Ft_22(I1)/(60*1000),fEPSP_22(I1)-fEPSP_keinDA_22(I1),'m^-',...
     Ft_23(I1)/(60*1000),fEPSP_23(I1)-fEPSP_keinDA_23(I1),'ro-',Ft_24(I1)/(60*1000),fEPSP_24(I1)-fEPSP_keinDA_24(I1),'b<-',...
     Ft_25(I1)/(60*1000),fEPSP_25(I1)-fEPSP_keinDA_25(I1),'g*-',Ft_26(I1)/(60*1000),fEPSP_26(I1)-fEPSP_keinDA_26(I1),'c>-',...
     Ft_27(I1)/(60*1000),fEPSP_27(I1)-fEPSP_keinDA_27(I1),'kd-',Ft_28(I1)/(60*1000),fEPSP_28(I1)-fEPSP_keinDA_28(I1),'m^-',...
     HFS,[55,62],'k-',SKF,[25,25],'b',...
     'LineWidth',ls,'MarkerSize',ms)
% plot(Ft_1/(60*1000)-2,fEPSP_1-fEPSP_keinDA_1,'ks-',FT_HFS_test_SKF50m_m200/(60*1000)-1,fEPSP_HFS_test_SKF50m_m200-fEPSP_HFS_test_keinDA_SKF50m_m200,'rs-','LineWidth',ls,'MarkerSize',ms)
% legend('100 \muM SKF 38393','50 \muM SKF 38393','25 \muM SKF 38393','15 \muM SKF 38393','10 \muM SKF 38393','5 \muM SKF 38393','2 \muM SKF 38393','1 \muM SKF 38393')
ylabel('\Delta fEPSP')
xlabel('Time (min)')
% title('\Deltat = -30')
xlim([0,225])
ylim([-20,65])
xticks([0:50:450])




HFS = 23*[1,1]; 
SKF = [3,23];

N1 = length(Ft_31);
a = 110;
b = 112;
I1 = [1:nopts:N1];
% I1 = cat(2,[1:nopts:a],[b:nopts:N1]);

% figure(3)%'DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
subplot(3,1,3)
plot(Ft_31(I1)/(60*1000),fEPSP_31(I1)-fEPSP_keinDA_31(I1),'ks-',...Ft_32(I1)/(60*1000),fEPSP_32(I1)-fEPSP_keinDA_32(I1),'m^-',...
     Ft_33(I1)/(60*1000),fEPSP_33(I1)-fEPSP_keinDA_33(I1),'ro-',Ft_34(I1)/(60*1000),fEPSP_34(I1)-fEPSP_keinDA_34(I1),'b<-',...
     Ft_35(I1)/(60*1000),fEPSP_35(I1)-fEPSP_keinDA_35(I1),'g*-',Ft_36(I1)/(60*1000),fEPSP_36(I1)-fEPSP_keinDA_36(I1),'c>-',...
     Ft_37(I1)/(60*1000),fEPSP_37(I1)-fEPSP_keinDA_37(I1),'kd-',Ft_38(I1)/(60*1000),fEPSP_38(I1)-fEPSP_keinDA_38(I1),'m^-',...
     HFS,[55,62],'k-',SKF,[25,25],'b',...
     'LineWidth',ls,'MarkerSize',ms)
% plot(Ft_1/(60*1000)-2,fEPSP_1-fEPSP_keinDA_1,'ks-',FT_HFS_test_SKF50m_m200/(60*1000)-1,fEPSP_HFS_test_SKF50m_m200-fEPSP_HFS_test_keinDA_SKF50m_m200,'rs-','LineWidth',ls,'MarkerSize',ms)
% legend('100 \muM SKF 38393','50 \muM SKF 38393','25 \muM SKF 38393','15 \muM SKF 38393','10 \muM SKF 38393','5 \muM SKF 38393','2 \muM SKF 38393','1 \muM SKF 38393')
ylabel('\Delta fEPSP')
xlabel('Time (min)')
% title('\Deltat = -20')
% ylim([-20,65])
xlim([0,225])
ylim([-20,65])
xticks([0:50:300])

