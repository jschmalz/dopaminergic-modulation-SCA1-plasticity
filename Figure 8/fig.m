clearvars
clc
clear -all

load('t_SKF50_m212.mat')
load('t_SKF50_m212_keinDA.mat')
load('fEPSP_SKF50_m212_keinDA.mat')
load('fEPSP_SKF50_m212.mat')

fEPSP_1 = fEPSP_SKF50_m212;
Ft_1 = t_SKF50_m212;
fEPSP_keinDA_1 = fEPSP_SKF50_m212_keinDA;
Ft_keinDA_1 = t_SKF50_m212_keinDA;          

load('t_SKF25_SHFS_m212.mat')
load('t_SKF25_SHFS_m212_keinDA.mat')
load('fEPSP_SKF25_SHFS_m212.mat')
load('fEPSP_SKF25_SHFS_m212_keinDA.mat')

fEPSP_2 = fEPSP_SKF25_SHFS_m212;
Ft_2 = t_SKF25_SHFS_m212;
fEPSP_keinDA_2 = fEPSP_SKF25_SHFS_m212_keinDA;
Ft_keinDA_2 = t_SKF25_SHFS_m212_keinDA;          

load('t_SKF15_SHFS_m212.mat')
load('t_SKF15_SHFS_m212_keinDA.mat')
load('fEPSP_SKF15_SHFS_m212.mat')
load('fEPSP_SKF15_SHFS_m212_keinDA.mat')

fEPSP_3 = fEPSP_SKF15_SHFS_m212;
Ft_3 = t_SKF15_SHFS_m212;
fEPSP_keinDA_3 = fEPSP_SKF15_SHFS_m212_keinDA;
Ft_keinDA_3 = t_SKF15_SHFS_m212_keinDA;          

load('t_SKF10_SHFS_m212.mat')
load('t_SKF10_SHFS_m212_keinDA.mat')
load('fEPSP_SKF10_SHFS_m212.mat')
load('fEPSP_SKF10_SHFS_m212_keinDA.mat')

fEPSP_4 = fEPSP_SKF10_SHFS_m212;
Ft_4 = t_SKF10_SHFS_m212;
fEPSP_keinDA_4 = fEPSP_SKF10_SHFS_m212_keinDA;
Ft_keinDA_4 = t_SKF10_SHFS_m212_keinDA;          

load('t_SKF5_SHFS_m212.mat')
load('t_SKF5_SHFS_m212_keinDA.mat')
load('fEPSP_SKF5_SHFS_m212.mat')
load('fEPSP_SKF5_SHFS_m212_keinDA.mat')

fEPSP_5 = fEPSP_SKF5_SHFS_m212;
Ft_5 = t_SKF5_SHFS_m212;
fEPSP_keinDA_5 = fEPSP_SKF5_SHFS_m212_keinDA;
Ft_keinDA_5 = t_SKF5_SHFS_m212_keinDA;          

load('t_SKF2_SHFS_m212.mat')
load('t_SKF2_SHFS_m212_keinDA.mat')
load('fEPSP_SKF2_SHFS_m212.mat')
load('fEPSP_SKF2_SHFS_m212_keinDA.mat')

fEPSP_6 = fEPSP_SKF2_SHFS_m212;
Ft_6 = t_SKF2_SHFS_m212;
fEPSP_keinDA_6 = fEPSP_SKF2_SHFS_m212_keinDA;
Ft_keinDA_6 = t_SKF2_SHFS_m212_keinDA;          

load('t_SKF1_SHFS_m212.mat')
load('t_SKF1_SHFS_m212_keinDA.mat')
load('fEPSP_SKF1_SHFS_m212.mat')
load('fEPSP_SKF1_SHFS_m212_keinDA.mat')

fEPSP_7 = fEPSP_SKF1_SHFS_m212;
Ft_7 = t_SKF1_SHFS_m212;
fEPSP_keinDA_7 = fEPSP_SKF1_SHFS_m212_keinDA;
Ft_keinDA_7 = t_SKF1_SHFS_m212_keinDA;          


% dt = -30

load('t_SKF50_SHFS_m30.mat')
load('t_SKF50_SHFS_m30_keinDA.mat')
load('fEPSP_SKF50_SHFS_m30.mat')
load('fEPSP_SKF50_SHFS_m30_keinDA.mat')

fEPSP_21 = fEPSP_SKF50_SHFS_m30;
Ft_21 = t_SKF50_SHFS_m30;
fEPSP_keinDA_21 = fEPSP_SKF50_SHFS_m30_keinDA;
Ft_keinDA_21 = t_SKF50_SHFS_m30_keinDA; 


load('t_SKF25_SHFS_m30.mat')
load('t_SKF25_SHFS_m30_keinDA.mat')
load('fEPSP_SKF25_SHFS_m30.mat')
load('fEPSP_SKF25_SHFS_m30_keinDA.mat')

fEPSP_22 = fEPSP_SKF25_SHFS_m30;
Ft_22 = t_SKF25_SHFS_m30;
fEPSP_keinDA_22 = fEPSP_SKF25_SHFS_m30_keinDA;
Ft_keinDA_22 = t_SKF25_SHFS_m30_keinDA; 


load('t_SKF15_SHFS_m30.mat')
load('t_SKF15_SHFS_m30_keinDA.mat')
load('fEPSP_SKF15_SHFS_m30.mat')
load('fEPSP_SKF15_SHFS_m30_keinDA.mat')

fEPSP_23 = fEPSP_SKF15_SHFS_m30;
Ft_23 = t_SKF15_SHFS_m30;
fEPSP_keinDA_23 = fEPSP_SKF15_SHFS_m30_keinDA;
Ft_keinDA_23 = t_SKF15_SHFS_m30_keinDA;


load('t_SKF10_SHFS_m30.mat')
load('t_SKF10_SHFS_m30_keinDA.mat')
load('fEPSP_SKF10_SHFS_m30.mat')
load('fEPSP_SKF10_SHFS_m30_keinDA.mat')

fEPSP_24 = fEPSP_SKF10_SHFS_m30;
Ft_24 = t_SKF10_SHFS_m30;
fEPSP_keinDA_24 = fEPSP_SKF10_SHFS_m30_keinDA;
Ft_keinDA_24 = t_SKF10_SHFS_m30_keinDA;

load('t_SKF5_SHFS_m30.mat')
load('t_SKF5_SHFS_m30_keinDA.mat')
load('fEPSP_SKF5_SHFS_m30.mat')
load('fEPSP_SKF5_SHFS_m30_keinDA.mat')
% 
fEPSP_25 = fEPSP_SKF5_SHFS_m30;
Ft_25 = t_SKF5_SHFS_m30;
fEPSP_keinDA_25 = fEPSP_SKF5_SHFS_m30_keinDA;
Ft_keinDA_25 = t_SKF5_SHFS_m30_keinDA;


load('t_SKF2_SHFS_m30.mat')
load('t_SKF2_SHFS_m30_keinDA.mat')
load('fEPSP_SKF2_SHFS_m30.mat')
load('fEPSP_SKF2_SHFS_m30_keinDA.mat')

fEPSP_26 = fEPSP_SKF2_SHFS_m30;
Ft_26 = t_SKF2_SHFS_m30;
fEPSP_keinDA_26 = fEPSP_SKF2_SHFS_m30_keinDA;
Ft_keinDA_26 = t_SKF2_SHFS_m30_keinDA;


load('t_SKF1_SHFS_m30.mat')
load('t_SKF1_SHFS_m30_keinDA.mat')
load('fEPSP_SKF1_SHFS_m30.mat')
load('fEPSP_SKF1_SHFS_m30_keinDA.mat')

fEPSP_27 = fEPSP_SKF1_SHFS_m30;
Ft_27 = t_SKF1_SHFS_m30;
fEPSP_keinDA_27 = fEPSP_SKF1_SHFS_m30_keinDA;
Ft_keinDA_27 = t_SKF1_SHFS_m30_keinDA;

% dt = -15

load('t_SKF50_SHFS_m15.mat')
load('t_SKF50_SHFS_m15_keinDA.mat')
load('fEPSP_SKF50_SHFS_m15.mat')
load('fEPSP_SKF50_SHFS_m15_keinDA.mat')

fEPSP_31 = fEPSP_SKF50_SHFS_m15;
Ft_31 = t_SKF50_SHFS_m15;
fEPSP_keinDA_31 = fEPSP_SKF50_SHFS_m15_keinDA;
Ft_keinDA_31 = t_SKF50_SHFS_m15_keinDA; 


load('t_SKF25_SHFS_m15.mat')
load('t_SKF25_SHFS_m15_keinDA.mat')
load('fEPSP_SKF25_SHFS_m15.mat')
load('fEPSP_SKF25_SHFS_m15_keinDA.mat')

fEPSP_32 = fEPSP_SKF25_SHFS_m15;
Ft_32 = t_SKF25_SHFS_m15;
fEPSP_keinDA_32 = fEPSP_SKF25_SHFS_m15_keinDA;
Ft_keinDA_32 = t_SKF25_SHFS_m15_keinDA; 


load('t_SKF15_SHFS_m15.mat')
load('t_SKF15_SHFS_m15_keinDA.mat')
load('fEPSP_SKF15_SHFS_m15.mat')
load('fEPSP_SKF15_SHFS_m15_keinDA.mat')

fEPSP_33 = fEPSP_SKF15_SHFS_m15;
Ft_33 = t_SKF15_SHFS_m15;
fEPSP_keinDA_33 = fEPSP_SKF15_SHFS_m15_keinDA;
Ft_keinDA_33 = t_SKF15_SHFS_m15_keinDA;


load('t_SKF10_SHFS_m15.mat')
load('t_SKF10_SHFS_m15_keinDA.mat')
load('fEPSP_SKF10_SHFS_m15.mat')
load('fEPSP_SKF10_SHFS_m15_keinDA.mat')

fEPSP_34 = fEPSP_SKF10_SHFS_m15;
Ft_34 = t_SKF10_SHFS_m15;
fEPSP_keinDA_34 = fEPSP_SKF10_SHFS_m15_keinDA;
Ft_keinDA_34 = t_SKF10_SHFS_m15_keinDA;

load('t_SKF5_SHFS_m15.mat')
load('t_SKF5_SHFS_m15_keinDA.mat')
load('fEPSP_SKF5_SHFS_m15.mat')
load('fEPSP_SKF5_SHFS_m15_keinDA.mat')
% 
fEPSP_35 = fEPSP_SKF5_SHFS_m15;
Ft_35 = t_SKF5_SHFS_m15;
fEPSP_keinDA_35 = fEPSP_SKF5_SHFS_m15_keinDA;
Ft_keinDA_35 = t_SKF5_SHFS_m15_keinDA;


load('t_SKF2_SHFS_m15.mat')
load('t_SKF2_SHFS_m15_keinDA.mat')
load('fEPSP_SKF2_SHFS_m15.mat')
load('fEPSP_SKF2_SHFS_m15_keinDA.mat')

fEPSP_36 = fEPSP_SKF2_SHFS_m15;
Ft_36 = t_SKF2_SHFS_m15;
fEPSP_keinDA_36 = fEPSP_SKF2_SHFS_m15_keinDA;
Ft_keinDA_36 = t_SKF2_SHFS_m15_keinDA;


load('t_SKF1_SHFS_m15.mat')
load('t_SKF1_SHFS_m15_keinDA.mat')
load('fEPSP_SKF1_SHFS_m15.mat')
load('fEPSP_SKF1_SHFS_m15_keinDA.mat')

fEPSP_37 = fEPSP_SKF1_SHFS_m15;
Ft_37 = t_SKF1_SHFS_m15;
fEPSP_keinDA_37 = fEPSP_SKF1_SHFS_m15_keinDA;
Ft_keinDA_37 = t_SKF1_SHFS_m15_keinDA;




ls = 1.5;
ms = 6;

HFS = 214*[1,1]; 
SKF = [2,17];

nopts = 5;

N1 = length(Ft_1);
I1 = [1:nopts:N1];

figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
subplot(3,1,1)
plot(Ft_1(I1)/(60*1000),fEPSP_1(I1)-fEPSP_keinDA_1(I1),'ks-',Ft_2(I1)/(60*1000),fEPSP_2(I1)-fEPSP_keinDA_2(I1),'m^-',Ft_3(I1)/(60*1000),fEPSP_3(I1)-fEPSP_keinDA_3(I1),'ro-',...
    Ft_4(I1)/(60*1000),fEPSP_4(I1)-fEPSP_keinDA_4(I1),'b<-',Ft_5(I1)/(60*1000),fEPSP_5(I1)-fEPSP_keinDA_5(I1),'g*-',...
    Ft_6(I1)/(60*1000),fEPSP_6(I1)-fEPSP_keinDA_6(I1),'c>-',Ft_7(I1)/(60*1000),fEPSP_7(I1)-fEPSP_keinDA_7(I1),'kd-',...
    HFS,[52,63],'k-',HFS+10,[52,63],'k-',HFS+20,[52,63],'k-',SKF,[25,25],'b',...
    'LineWidth',ls,'MarkerSize',ms)
% plot(Ft_1/(60*1000)-2,fEPSP_1-fEPSP_keinDA_1,'ks-',FT_HFS_test_SKF50m_m200/(60*1000)-1,fEPSP_HFS_test_SKF50m_m200-fEPSP_HFS_test_keinDA_SKF50m_m200,'rs-','LineWidth',ls,'MarkerSize',ms)
legend('50 \muM SKF 38393','25 \muM SKF 38393','15 \muM SKF 38393','10 \muM SKF 38393','5 \muM SKF 38393','2 \muM SKF 38393','1 \muM SKF 38393')
ylabel('\Delta fEPSP')
xlabel('Time (min)')
title('\Deltat = -212')
xlim([0,320])
ylim([-15,65])
xticks([0:50:300])

N2 = length(Ft_21);
I2 = [1:nopts:N2];
HFS = 32*[1,1]; 
SKF = [2,17];



% figure(2)
subplot(3,1,2)
plot(Ft_21(I2)/(60*1000),fEPSP_21(I2)-fEPSP_keinDA_21(I2),'ks-',Ft_22(I2)/(60*1000),fEPSP_22(I2)-fEPSP_keinDA_22(I2),'m^-',Ft_23(I2)/(60*1000),fEPSP_23(I2)-fEPSP_keinDA_23(I2),'ro-',...
     Ft_24(I2)/(60*1000),fEPSP_24(I2)-fEPSP_keinDA_24(I2),'b<-',Ft_25(I2)/(60*1000),fEPSP_25(I2)-fEPSP_keinDA_25(I2),'g*-',...
    Ft_26(I2)/(60*1000),fEPSP_26(I2)-fEPSP_keinDA_26(I2),'c>-',Ft_27(I2)/(60*1000),fEPSP_27(I2)-fEPSP_keinDA_27(I2),'kd-',...
    HFS,[30,39],'k-',HFS+10,[30,39],'k-',HFS+20,[30,39],'k-',SKF,[25,25],'b',...
    'LineWidth',ls,'MarkerSize',ms)
% plot(Ft_1/(60*1000)-2,fEPSP_1-fEPSP_keinDA_1,'ks-',FT_HFS_test_SKF50m_m200/(60*1000)-1,fEPSP_HFS_test_SKF50m_m200-fEPSP_HFS_test_keinDA_SKF50m_m200,'rs-','LineWidth',ls,'MarkerSize',ms)
% legend('50 \muM SKF 38393','25 \muM SKF 38393','15 \muM SKF 38393','10 \muM SKF 38393','5 \muM SKF 38393','2 \muM SKF 38393','1 \muM SKF 38393')
ylabel('\Delta fEPSP')
xlabel('Time (min)')
title('\Deltat = -30')
xlim([0,220])
ylim([-15,40])
xticks([0:50:200])


N3 = length(Ft_31);
I3 = [1:nopts:N3];

HFS = 17*[1,1]; 
SKF = [2,17];


% figure(3)
subplot(3,1,3)
plot(Ft_31(I3)/(60*1000),fEPSP_31(I3)-fEPSP_keinDA_31(I3),'ks-',Ft_32(I3)/(60*1000),fEPSP_32(I3)-fEPSP_keinDA_32(I3),'m^-',Ft_33(I3)/(60*1000),fEPSP_33(I3)-fEPSP_keinDA_33(I3),'ro-',...
     Ft_34(I3)/(60*1000),fEPSP_34(I3)-fEPSP_keinDA_34(I3),'b<-',Ft_35(I3)/(60*1000),fEPSP_35(I3)-fEPSP_keinDA_35(I3),'g*-',...
    Ft_36(I3)/(60*1000),fEPSP_36(I3)-fEPSP_keinDA_36(I3),'c>-',Ft_37(I3)/(60*1000),fEPSP_37(I3)-fEPSP_keinDA_37(I3),'kd-',...
    HFS,[30,39],'k-',HFS+10,[30,39],'k-',HFS+20,[30,39],'k-',SKF,[25,25],'b',...
    'LineWidth',ls,'MarkerSize',ms)
% plot(Ft_1/(60*1000)-2,fEPSP_1-fEPSP_keinDA_1,'ks-',FT_HFS_test_SKF50m_m200/(60*1000)-1,fEPSP_HFS_test_SKF50m_m200-fEPSP_HFS_test_keinDA_SKF50m_m200,'rs-','LineWidth',ls,'MarkerSize',ms)
% legend('50 \muM SKF 38393','25 \muM SKF 38393','15 \muM SKF 38393','10 \muM SKF 38393','5 \muM SKF 38393','2 \muM SKF 38393','1 \muM SKF 38393')
ylabel('\Delta fEPSP')
xlabel('Time (min)')
title('\Deltat = -15')
xlim([0,210])
ylim([-15,40])
xticks([0:50:200])

