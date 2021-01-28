clearvars
clear -all
clc

% dt = 10


load('t_SKF50_WHFS_10.mat')
load('t_SKF50_WHFS_10_keinDA.mat')
load('fEPSP_SKF50_WHFS_10.mat')
load('fEPSP_SKF50_WHFS_10_keinDA.mat')

fEPSP_41 = fEPSP_SKF50_WHFS_10;
Ft_41 = t_SKF50_WHFS_10;
fEPSP_keinDA_41 = fEPSP_SKF50_WHFS_10_keinDA;
Ft_keinDA_41 = t_SKF50_WHFS_10_keinDA; 


load('t_SKF25_WHFS_10.mat')
load('t_SKF25_WHFS_10_keinDA.mat')
load('fEPSP_SKF25_WHFS_10.mat')
load('fEPSP_SKF25_WHFS_10_keinDA.mat')

fEPSP_42 = fEPSP_SKF25_WHFS_10;
Ft_42 = t_SKF25_WHFS_10;
fEPSP_keinDA_42 = fEPSP_SKF25_WHFS_10_keinDA;
Ft_keinDA_42 = t_SKF25_WHFS_10_keinDA; 


load('t_SKF15_WHFS_10.mat')
load('t_SKF15_WHFS_10_keinDA.mat')
load('fEPSP_SKF15_WHFS_10.mat')
load('fEPSP_SKF15_WHFS_10_keinDA.mat')

fEPSP_43 = fEPSP_SKF15_WHFS_10;
Ft_43 = t_SKF15_WHFS_10;
fEPSP_keinDA_43 = fEPSP_SKF15_WHFS_10_keinDA;
Ft_keinDA_43 = t_SKF15_WHFS_10_keinDA;


load('t_SKF10_WHFS_10.mat')
load('t_SKF10_WHFS_10_keinDA.mat')
load('fEPSP_SKF10_WHFS_10.mat')
load('fEPSP_SKF10_WHFS_10_keinDA.mat')

fEPSP_44 = fEPSP_SKF10_WHFS_10;
Ft_44 = t_SKF10_WHFS_10;
fEPSP_keinDA_44 = fEPSP_SKF10_WHFS_10_keinDA;
Ft_keinDA_44 = t_SKF10_WHFS_10_keinDA;

load('t_SKF5_WHFS_10.mat')
load('t_SKF5_WHFS_10_keinDA.mat')
load('fEPSP_SKF5_WHFS_10.mat')
load('fEPSP_SKF5_WHFS_10_keinDA.mat')
% 
fEPSP_45 = fEPSP_SKF5_WHFS_10;
Ft_45 = t_SKF5_WHFS_10;
fEPSP_keinDA_45 = fEPSP_SKF5_WHFS_10_keinDA;
Ft_keinDA_45 = t_SKF5_WHFS_10_keinDA;


load('t_SKF2_WHFS_10.mat')
load('t_SKF2_WHFS_10_keinDA.mat')
load('fEPSP_SKF2_WHFS_10.mat')
load('fEPSP_SKF2_WHFS_10_keinDA.mat')

fEPSP_46 = fEPSP_SKF2_WHFS_10;
Ft_46 = t_SKF2_WHFS_10;
fEPSP_keinDA_46 = fEPSP_SKF2_WHFS_10_keinDA;
Ft_keinDA_46 = t_SKF2_WHFS_10_keinDA;


load('t_SKF1_WHFS_10.mat')
load('t_SKF1_WHFS_10_keinDA.mat')
load('fEPSP_SKF1_WHFS_10.mat')
load('fEPSP_SKF1_WHFS_10_keinDA.mat')

fEPSP_47 = fEPSP_SKF1_WHFS_10;
Ft_47 = t_SKF1_WHFS_10;
fEPSP_keinDA_47 = fEPSP_SKF1_WHFS_10_keinDA;
Ft_keinDA_47 = t_SKF1_WHFS_10_keinDA;


% dt = 30



load('t_SKF50_WHFS_30.mat')
load('t_SKF50_WHFS_30_keinDA.mat')
load('fEPSP_SKF50_WHFS_30.mat')
load('fEPSP_SKF50_WHFS_30_keinDA.mat')

fEPSP_51 = fEPSP_SKF50_WHFS_30;
Ft_51 = t_SKF50_WHFS_30;
fEPSP_keinDA_51 = fEPSP_SKF50_WHFS_30_keinDA;
Ft_keinDA_51 = t_SKF50_WHFS_30_keinDA; 


load('t_SKF25_WHFS_30.mat')
load('t_SKF25_WHFS_30_keinDA.mat')
load('fEPSP_SKF25_WHFS_30.mat')
load('fEPSP_SKF25_WHFS_30_keinDA.mat')

fEPSP_52 = fEPSP_SKF25_WHFS_30;
Ft_52 = t_SKF25_WHFS_30;
fEPSP_keinDA_52 = fEPSP_SKF25_WHFS_30_keinDA;
Ft_keinDA_52 = t_SKF25_WHFS_30_keinDA; 


load('t_SKF15_WHFS_30.mat')
load('t_SKF15_WHFS_30_keinDA.mat')
load('fEPSP_SKF15_WHFS_30.mat')
load('fEPSP_SKF15_WHFS_30_keinDA.mat')

fEPSP_53 = fEPSP_SKF15_WHFS_30;
Ft_53 = t_SKF15_WHFS_30;
fEPSP_keinDA_53 = fEPSP_SKF15_WHFS_30_keinDA;
Ft_keinDA_53 = t_SKF15_WHFS_30_keinDA;


load('t_SKF10_WHFS_30.mat')
load('t_SKF10_WHFS_30_keinDA.mat')
load('fEPSP_SKF10_WHFS_30.mat')
load('fEPSP_SKF10_WHFS_30_keinDA.mat')

fEPSP_54 = fEPSP_SKF10_WHFS_30;
Ft_54 = t_SKF10_WHFS_30;
fEPSP_keinDA_54 = fEPSP_SKF10_WHFS_30_keinDA;
Ft_keinDA_54 = t_SKF10_WHFS_30_keinDA;

load('t_SKF5_WHFS_30.mat')
load('t_SKF5_WHFS_30_keinDA.mat')
load('fEPSP_SKF5_WHFS_30.mat')
load('fEPSP_SKF5_WHFS_30_keinDA.mat')
% 
fEPSP_55 = fEPSP_SKF5_WHFS_30;
Ft_55 = t_SKF5_WHFS_30;
fEPSP_keinDA_55 = fEPSP_SKF5_WHFS_30_keinDA;
Ft_keinDA_55 = t_SKF5_WHFS_30_keinDA;


load('t_SKF2_WHFS_30.mat')
load('t_SKF2_WHFS_30_keinDA.mat')
load('fEPSP_SKF2_WHFS_30.mat')
load('fEPSP_SKF2_WHFS_30_keinDA.mat')

fEPSP_56 = fEPSP_SKF2_WHFS_30;
Ft_56 = t_SKF2_WHFS_30;
fEPSP_keinDA_56 = fEPSP_SKF2_WHFS_30_keinDA;
Ft_keinDA_56 = t_SKF2_WHFS_30_keinDA;


load('t_SKF1_WHFS_30.mat')
load('t_SKF1_WHFS_30_keinDA.mat')
load('fEPSP_SKF1_WHFS_30.mat')
load('fEPSP_SKF1_WHFS_30_keinDA.mat')

fEPSP_57 = fEPSP_SKF1_WHFS_30;
Ft_57 = t_SKF1_WHFS_30;
fEPSP_keinDA_57 = fEPSP_SKF1_WHFS_30_keinDA;
Ft_keinDA_57 = t_SKF1_WHFS_30_keinDA;

%%%%%%%%%%%%%%
% dt = 60
%%%%%%%%%%%%%%

load('t_SKF50_WHFS_60.mat')
load('t_SKF50_WHFS_60_keinDA.mat')
load('fEPSP_SKF50_WHFS_60.mat')
load('fEPSP_SKF50_WHFS_60keinDA.mat')

fEPSP_61 = fEPSP_SKF50_WHFS_60;
Ft_61 = t_SKF50_WHFS_60;
fEPSP_keinDA_61 = fEPSP_SKF50_WHFS_60_keinDA;
Ft_keinDA_61 = t_SKF50_WHFS_60_keinDA; 


load('t_SKF25_WHFS_60.mat')
load('t_SKF25_WHFS_60_keinDA.mat')
load('fEPSP_SKF25_WHFS_60.mat')
load('fEPSP_SKF25_WHFS_60keinDA.mat')

fEPSP_62 = fEPSP_SKF25_WHFS_60;
Ft_62 = t_SKF25_WHFS_60;
fEPSP_keinDA_62 = fEPSP_SKF25_WHFS_60_keinDA;
Ft_keinDA_62 = t_SKF25_WHFS_60_keinDA; 


load('t_SKF15_WHFS_60.mat')
load('t_SKF15_WHFS_60_keinDA.mat')
load('fEPSP_SKF15_WHFS_60.mat')
load('fEPSP_SKF15_WHFS_60keinDA.mat')

fEPSP_63 = fEPSP_SKF15_WHFS_60;
Ft_63 = t_SKF15_WHFS_60;
fEPSP_keinDA_63 = fEPSP_SKF15_WHFS_60_keinDA;
Ft_keinDA_63 = t_SKF15_WHFS_60_keinDA;


load('t_SKF10_WHFS_60.mat')
load('t_SKF10_WHFS_60_keinDA.mat')
load('fEPSP_SKF10_WHFS_60.mat')
load('fEPSP_SKF10_WHFS_60keinDA.mat')

fEPSP_64 = fEPSP_SKF10_WHFS_60;
Ft_64 = t_SKF10_WHFS_60;
fEPSP_keinDA_64 = fEPSP_SKF10_WHFS_60_keinDA;
Ft_keinDA_64 = t_SKF10_WHFS_60_keinDA;

load('t_SKF5_WHFS_60.mat')
load('t_SKF5_WHFS_60_keinDA.mat')
load('fEPSP_SKF5_WHFS_60.mat')
load('fEPSP_SKF5_WHFS_60keinDA.mat')
% 
fEPSP_65 = fEPSP_SKF5_WHFS_60;
Ft_65 = t_SKF5_WHFS_60;
fEPSP_keinDA_65 = fEPSP_SKF5_WHFS_60_keinDA;
Ft_keinDA_65 = t_SKF5_WHFS_60_keinDA;


load('t_SKF2_WHFS_60.mat')
load('t_SKF2_WHFS_60_keinDA.mat')
load('fEPSP_SKF2_WHFS_60.mat')
load('fEPSP_SKF2_WHFS_60keinDA.mat')

fEPSP_66 = fEPSP_SKF2_WHFS_60;
Ft_66 = t_SKF2_WHFS_60;
fEPSP_keinDA_66 = fEPSP_SKF2_WHFS_60_keinDA;
Ft_keinDA_66 = t_SKF2_WHFS_60_keinDA;


load('t_SKF1_WHFS_60.mat')
load('t_SKF1_WHFS_60_keinDA.mat')
load('fEPSP_SKF1_WHFS_60.mat')
load('fEPSP_SKF1_WHFS_60keinDA.mat')

fEPSP_67 = fEPSP_SKF1_WHFS_60;
Ft_67 = t_SKF1_WHFS_60;
fEPSP_keinDA_67 = fEPSP_SKF1_WHFS_60_keinDA;
Ft_keinDA_67 = t_SKF1_WHFS_60_keinDA;


nopts = 5;
ls = 1.5;
ms = 6;

HFS = 12*[1,1]; 
SKF = [(12+10),(12+10+15)];
N4 = length(Ft_41);
I4 = [1:nopts:N4];

% figure(10)
% plot(Ft_41/(60*1000),fEPSP_41,Ft_41/(60*1000),fEPSP_keinDA_41)

figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
% figure(1)
subplot(3,1,1)
plot(Ft_41(I4)/(60*1000),fEPSP_41(I4)-fEPSP_keinDA_41(I4),'ks-',Ft_42(I4)/(60*1000),fEPSP_42(I4)-fEPSP_keinDA_42(I4),'m^-',Ft_43(I4)/(60*1000),fEPSP_43(I4)-fEPSP_keinDA_43(I4),'ro-',...
     Ft_44(I4)/(60*1000),fEPSP_44(I4)-fEPSP_keinDA_44(I4),'b<-',Ft_45(I4)/(60*1000),fEPSP_45(I4)-fEPSP_keinDA_45(I4),'g*-',...
    Ft_46(I4)/(60*1000),fEPSP_46(I4)-fEPSP_keinDA_46(I4),'c>-',Ft_47(I4)/(60*1000),fEPSP_47(I4)-fEPSP_keinDA_47(I4),'kd-',...
    HFS,[5,8],'k-',SKF,[12,12],'b',...
    'LineWidth',ls,'MarkerSize',ms)
% plot(Ft_1/(60*1000)-2,fEPSP_1-fEPSP_keinDA_1,'ks-',FT_HFS_test_SKF50m_m200/(60*1000)-1,fEPSP_HFS_test_SKF50m_m200-fEPSP_HFS_test_keinDA_SKF50m_m200,'rs-','LineWidth',ls,'MarkerSize',ms)
legend('50 \muM SKF 38393','25 \muM SKF 38393','15 \muM SKF 38393','10 \muM SKF 38393','5 \muM SKF 38393','2 \muM SKF 38393','1 \muM SKF 38393')
ylabel('\Delta fEPSP')
xlabel('Time (min)')
title('\Deltat = 15')
ylim([-2,18])
xlim([0,220])
xticks([0:50:200])


N5 = length(Ft_51);
I5 = [1:nopts:N5];

HFS = 12*[1,1]; 
SKF = [(12+30),(12+30+15)];

subplot(3,1,2)
% figure(2)
plot(Ft_51(I5)/(60*1000),fEPSP_51(I5)-fEPSP_keinDA_51(I5),'ks-',Ft_52(I5)/(60*1000),fEPSP_52(I5)-fEPSP_keinDA_52(I5),'m^-',Ft_53(I5)/(60*1000),fEPSP_53(I5)-fEPSP_keinDA_53(I5),'ro-',...
     Ft_54(I5)/(60*1000),fEPSP_54(I5)-fEPSP_keinDA_54(I5),'b<-',Ft_55(I5)/(60*1000),fEPSP_55(I5)-fEPSP_keinDA_55(I5),'g*-',...
    Ft_56(I5)/(60*1000),fEPSP_56(I5)-fEPSP_keinDA_56(I5),'c>-',Ft_57(I5)/(60*1000),fEPSP_57(I5)-fEPSP_keinDA_57(I5),'kd-',...
     HFS,[5,8],'k-',SKF,[12,12],'b',...
    'LineWidth',ls,'MarkerSize',ms)
% plot(Ft_1/(60*1000)-2,fEPSP_1-fEPSP_keinDA_1,'ks-',FT_HFS_test_SKF50m_m200/(60*1000)-1,fEPSP_HFS_test_SKF50m_m200-fEPSP_HFS_test_keinDA_SKF50m_m200,'rs-','LineWidth',ls,'MarkerSize',ms)
% legend('50 \muM SKF 38393','25 \muM SKF 38393','15 \muM SKF 38393','10 \muM SKF 38393','5 \muM SKF 38393','2 \muM SKF 38393','1 \muM SKF 38393')
ylabel('\Delta fEPSP')
xlabel('Time (min)')
title('\Deltat = 30')
ylim([-2,18])
xlim([0,220])
xticks([0:50:200])

N6 = length(Ft_61);
I6 = [1:nopts:N6];

HFS = 12*[1,1]; 
SKF = [(12+60),(12+60+15)];

% figure(3)
subplot(3,1,3)
plot(Ft_61(I6)/(60*1000),fEPSP_61(I6)-fEPSP_keinDA_61(I6),'ks-',Ft_62(I6)/(60*1000),fEPSP_62(I6)-fEPSP_keinDA_62(I6),'m^-',Ft_63(I6)/(60*1000),fEPSP_63(I6)-fEPSP_keinDA_63(I6),'ro-',...
     Ft_64(I6)/(60*1000),fEPSP_64(I6)-fEPSP_keinDA_64(I6),'b<-',Ft_65(I6)/(60*1000),fEPSP_65(I6)-fEPSP_keinDA_65(I6),'g*-',...
    Ft_66(I6)/(60*1000),fEPSP_66(I6)-fEPSP_keinDA_66(I6),'c>-',Ft_67(I6)/(60*1000),fEPSP_67(I6)-fEPSP_keinDA_67(I6),'kd-',...
    HFS,[5,8],'k-',SKF,[12,12],'b',...
    'LineWidth',ls,'MarkerSize',ms)
% plot(Ft_1/(60*1000)-2,fEPSP_1-fEPSP_keinDA_1,'ks-',FT_HFS_test_SKF50m_m200/(60*1000)-1,fEPSP_HFS_test_SKF50m_m200-fEPSP_HFS_test_keinDA_SKF50m_m200,'rs-','LineWidth',ls,'MarkerSize',ms)
% legend('50 \muM SKF 38393','25 \muM SKF 38393','15 \muM SKF 38393','10 \muM SKF 38393','5 \muM SKF 38393','2 \muM SKF 38393','1 \muM SKF 38393')
ylabel('\Delta fEPSP')
xlabel('Time (min)')
title('\Deltat = 60')
ylim([-2,18])
xlim([0,270])
xticks([0:50:200])


