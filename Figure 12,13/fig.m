clearvars
clc
clear -all


load('FT_LFS3Hz_1200p_SKF100m_0.mat')
load('fEPSP_LFS3Hz_1200p_SKF100m_0.mat')
load('fEPSP_LFS3Hz_1200p_SKF100m_0_keinDA.mat')
load('FT_LFS3Hz_1200p_SKF100m_0_keinDA.mat')


t_1 = FT_LFS3Hz_1200p_SKF100m_0;
fEPSP_1 = fEPSP_LFS3Hz_1200p_SKF100m_0;
t_1_keinDA = FT_LFS3Hz_1200p_SKF100m_0_keinDA;
fEPSP_1_keinDA = fEPSP_LFS3Hz_1200p_SKF100m_0_keinDA;


load('FT_LFS3Hz_1200p_SKF100m_60.mat')
load('fEPSP_LFS3Hz_1200p_SKF100m_60.mat')
load('fEPSP_LFS3Hz_1200p_SKF100m_60_keinDA.mat')
load('FT_LFS3Hz_1200p_SKF100m_60_keinDA.mat')


t_2 = FT_LFS3Hz_1200p_SKF100m_60;
fEPSP_2 = fEPSP_LFS3Hz_1200p_SKF100m_60;
t_2_keinDA = FT_LFS3Hz_1200p_SKF100m_60_keinDA;
fEPSP_2_keinDA = fEPSP_LFS3Hz_1200p_SKF100m_60_keinDA;


load('fEPSP_LFS3Hz_1200p_SKF100m_m100.mat')
load('fEPSP_LFS3Hz_1200p_SKF100m_m100_keinDA.mat')
load('FT_LFS3Hz_1200p_SKF100m_m100.mat')
load('FT_LFS3Hz_1200p_SKF100m_m100_keinDA.mat')

t_3 = FT_LFS3Hz_1200p_SKF100m_m100;
fEPSP_3 = fEPSP_LFS3Hz_1200p_SKF100m_m100;
t_3_keinDA = FT_LFS3Hz_1200p_SKF100m_m100_keinDA;
fEPSP_3_keinDA = fEPSP_LFS3Hz_1200p_SKF100m_m100_keinDA;


load('fEPSP_LFS3Hz_1200p_SKF100m_m30.mat')
load('fEPSP_LFS3Hz_1200p_SKF100m_m30_keinDA.mat')
load('FT_LFS3Hz_1200p_SKF100m_m30.mat')
load('FT_LFS3Hz_1200p_SKF100m_m30_keinDA.mat')


t_4 = FT_LFS3Hz_1200p_SKF100m_m30;
fEPSP_4 = fEPSP_LFS3Hz_1200p_SKF100m_m30;
t_4_keinDA = FT_LFS3Hz_1200p_SKF100m_m30_keinDA;
fEPSP_4_keinDA = fEPSP_LFS3Hz_1200p_SKF100m_m30_keinDA;


load('mockett_2006_data.mat')
load('mockett_2007_2_data.mat')
load('Mockett_2007_LTD_lateSKF.mat')
load('Mockett_2007_LTD.mat')

ls = 1.5;
ms = 6;

HFS = 33*[1,1]; 
SKF = [0,20]+33+6;

nopts = 5;

N1 = length(t_1);


%% Figure 12
a = 1;
b = 17;
c = 21;



figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
subplot(2,1,1)
plot(Mockett_2007_LTD(:,1)+3,Mockett_2007_LTD(:,2)+100,'go',mockett_2006_data(:,1)+3,mockett_2006_data(:,2)+100,'ro',...
    t_1_keinDA(a:b)./60000,fEPSP_1_keinDA(a:b),'ks-',t_1(a:b)./60000,fEPSP_1(a:b),'b^-',...
    HFS,[110,120],'k-',SKF,[115,115],'b',...
    t_1_keinDA(c:end)./60000,fEPSP_1_keinDA(c:end),'ks-',t_1(c:end)./60000,fEPSP_1(c:end),'b^-',...
    'LineWidth',ls,'MarkerSize',ms)
ylabel('fEPSP (% change)')
xlabel('Time (min)')
xlim([0,160])
ylim([50,130])
legend('(Mockett, 2007) - LFS','(Mockett, 2007) - LFS + SKF','(Model) - LFS','(Model) - LFS + SKF')


%% Figure 13

a = 1;
b = 17;
c = 21;

HFS = 33*[1,1]; 
SKF = [0,20]+33+6+60;

subplot(2,1,2)
plot(Mockett_2007_LTD(:,1)+3,Mockett_2007_LTD(:,2)+100,'go',mockett_2007_2_data(:,1)+3,mockett_2007_2_data(:,2)+100,'ro',...
    t_2_keinDA(a:b)./60000,fEPSP_2_keinDA(a:b),'ks-',t_2(a:b)./60000,fEPSP_2(a:b),'b^-',...
    HFS,[110,120],'k-',SKF,[115,115],'b',...
    t_2_keinDA(c:end)./60000,fEPSP_2_keinDA(c:end),'ks-',t_2(c:end)./60000,fEPSP_2(c:end),'b^-',...
    'LineWidth',ls,'MarkerSize',ms)
ylabel('fEPSP (% change)')
xlabel('Time (min)')
xlim([0,160])
ylim([50,130])
legend('(Mockett, 2007) - LFS','(Mockett, 2007) - LFS + SKF','(Model) - LFS','(Model) - LFS + SKF')


%%

a = 1;
b = 52;
c = 56;

HFS = 103*[1,1]; 
SKF = [0,20]+3;



figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
subplot(2,1,1)
plot(t_3_keinDA(a:b)./60000,fEPSP_3_keinDA(a:b),'ks-',t_3(a:b)./60000,fEPSP_3(a:b),'b^-',...
    HFS,[155,170],'k-',SKF,[125,125],'b',...
    t_3_keinDA(c:end)./60000,fEPSP_3_keinDA(c:end),'ks-',t_3(c:end)./60000,fEPSP_3(c:end),'b^-',...
    'LineWidth',ls,'MarkerSize',ms)
ylabel('fEPSP (% change)')
xlabel('Time (min)')
xlim([0,200])
ylim([45,175])
legend('(Model) - LFS','(Model) - LFS + SKF')


a = 1;
b = 16;
c = 21;

HFS = 33*[1,1]; 
SKF = [0,20]+3;

subplot(2,1,2)
plot(t_4_keinDA(a:b)./60000,fEPSP_4_keinDA(a:b),'ks-',t_4(a:b)./60000,fEPSP_4(a:b),'b^-',...
    HFS,[135,145],'k-',SKF,[125,125],'b',...
    t_4_keinDA(c:end)./60000,fEPSP_4_keinDA(c:end),'ks-',t_4(c:end)./60000,fEPSP_4(c:end),'b^-',...
    'LineWidth',ls,'MarkerSize',ms)
ylabel('fEPSP (% change)')
xlabel('Time (min)')
xlim([0,200])
ylim([45,175])
legend('(Model) - LFS','(Model) - LFS + SKF')
