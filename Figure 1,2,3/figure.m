clearvars
clc
clear -all



% data
load('Huang_SKF_before.mat')
load('Huang_SKF_after.mat')
load('Navakkode_2012.mat')


td_1 = Huang_SKF_before(:,1);
xd_1 = Huang_SKF_before(:,2);
xd_1(14:end) = Huang_SKF_before(14:end,2) + Huang_SKF_before(13,2)-100;

td_2 = Huang_SKF_after(:,1);
xd_2 = Huang_SKF_after(:,2);
xd_2(9:end) = Huang_SKF_after(9:end,2) + mean(Huang_SKF_after(7:8,2)-100);

td_3 = Navakkode_2012(:,1);
xd_3 = Navakkode_2012(:,2);
xd_3(19:end) = Navakkode_2012(19:end,2) + Navakkode_2012(18,2)-100;



% model
load('t_SKF50_m212.mat')
load('t_SKF50_m212_keinDA.mat')
load('fEPSP_SKF50_m212_keinDA.mat')
load('fEPSP_SKF50_m212.mat')

fEPSP = fEPSP_SKF50_m212;
Ft = t_SKF50_m212;
fEPSP_keinDA = fEPSP_SKF50_m212_keinDA;
Ft_keinDA = t_SKF50_m212_keinDA;

load('t_SKF50_70.mat')
load('t_SKF50_70_keinDA.mat')
load('fEPSP_SKF50_70_keinDA.mat')
load('fEPSP_SKF50_70.mat')

fEPSP_2 = fEPSP_SKF50_70;
Ft_2 = t_SKF50_70;
fEPSP_keinDA_2 = fEPSP_SKF50_70_keinDA;
Ft_keinDA_2 = t_SKF50_70_keinDA;



load('fEPSP_DA50_m165.mat')
load('fEPSP_DA50_m165_keinDA.mat')
load('t_DA50_m165.mat')
load('t_DA50_m165_keinDA.mat')

fEPSP_3 = fEPSP_DA50_m165;
Ft_3 = t_DA50_m165;
fEPSP_keinDA_3 = fEPSP_DA50_m165_keinDA;
Ft_keinDA_3 = t_DA50_m165_keinDA;


load('t_SKF50_m30.mat')
load('t_SKF50_m30_keinDA.mat')
load('fEPSP_SKF50_m30.mat')
load('fEPSP_SKF50_m30_keinDA.mat')

fEPSP_4 = fEPSP_SKF50_m30;
Ft_4 = t_SKF50_m30;
fEPSP_keinDA_4 = fEPSP_SKF50_m30_keinDA;
Ft_keinDA_4 = t_SKF50_m30_keinDA;

load('t_SKF50_30.mat')
load('t_SKF50_30_keinDA.mat')
load('fEPSP_SKF50_30.mat')
load('fEPSP_SKF50_30_keinDA.mat')

fEPSP_5 = fEPSP_SKF50_30;
Ft_5 = t_SKF50_30;
fEPSP_keinDA_5 = fEPSP_SKF50_30_keinDA;
Ft_keinDA_5 = t_SKF50_30_keinDA;


ls = 1.5;
ms = 5;
D1RA = 0*[1,1];
STIM = 212*[1,1];


%% Figure 1

figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')%1)
subplot(2,1,1)
plot(Ft_keinDA./60000-2, fEPSP_keinDA,'ks-',Ft./60000-2,fEPSP,'^b-',td_1,xd_1,'ro',...
    D1RA,[140,170],'b',STIM,[270,300],'k',STIM+10,[270,300],'k',STIM+20,[270,300],'k','LineWidth',ls,'MarkerSize',ms)
xlabel('Time (sec)')
ylabel('fEPSP (% change)')
legend('(Model) - HFS','(Model) - SKF + HFS','(Huang, 1995) - SKF + HFS')
xlim([-2,400])
ylim([90,310])


D1RA = 70*[1,1];
STIM = 0*[1,1];

subplot(2,1,2)
plot(Ft_keinDA_2./60000-22, fEPSP_keinDA_2,'ks-',Ft_2./60000-22, fEPSP_2,'^b-',td_2,xd_2,'ro',...
    D1RA,[230,260],'b',STIM,[270,300],'k',STIM+10,[270,300],'k',STIM+20,[270,300],'k','LineWidth',ls,'MarkerSize',ms)
legend('(Model) - HFS','(Model) - HFS + SKF','(Huang, 1995) - HFS + SKF')
xlabel('Time (sec)')
ylabel('fEPSP (% change)')
xlim([-25,375])
ylim([90,310])

%% Figure 2

D1RA = 0*[1,1];
STIM = 165*[1,1];

figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')%1)
plot(Ft_keinDA_3./60000-2, fEPSP_keinDA_3,'ks-',Ft_3./60000-2, fEPSP_3,'^b-',td_3,xd_3,'ro',...
    D1RA,[140,170],'b',D1RA+10,[140,170],'b',D1RA+20,[140,170],'b',STIM,[270,300],'k',STIM+5,[270,300],'k',STIM+10,[270,300],'k','LineWidth',ls,'MarkerSize',ms)
legend('(Model) - HFS','(Model) - DA + HFS','(Navakkode, 2012) - DA + HFS')
xlabel('Time (sec)')
ylabel('fEPSP (% change)')
xlim([-2,400])
ylim([90,310])



%% Figure 3


D1RA = 0*[1,1]+2;
STIM = 30*[1,1]+2;


figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')%2)
subplot(2,1,1)
plot(Ft_keinDA_4./60000, fEPSP_keinDA_4,'ks-',Ft_4./60000, fEPSP_4,'^b-',...
    D1RA,[140,170],'b',STIM,[270,300],'k',STIM+10,[270,300],'k',STIM+20,[270,300],'k','LineWidth',ls,'MarkerSize',ms)
legend('(Model) - HFS','(Model) - SKF + HFS')
xlabel('Time (sec)')
ylabel('fEPSP (% change)')
ylim([90,310])
xlim([-2,400])

D1RA = 30*[1,1]+22;
STIM = 0*[1,1]+22;

subplot(2,1,2)
plot(Ft_keinDA_5./60000, fEPSP_keinDA_5,'ks-',Ft_5./60000, fEPSP_5,'^b-',...
    D1RA,[220,250],'b',STIM,[270,300],'k',STIM+10,[270,300],'k',STIM+20,[270,300],'k','LineWidth',ls,'MarkerSize',ms)
legend('(Model) - HFS','(Model) - HFS + SKF')
xlabel('Time (sec)')
ylabel('fEPSP (% change)')
ylim([90,310])
xlim([-2,400])





