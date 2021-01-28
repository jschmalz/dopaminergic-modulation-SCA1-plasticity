clearvars
clc
clear -all

load('FT_LFS1Hz_450p_SKF3m_m5.mat')
load('FT_LFS1Hz_450p_SKF3m_m5_keinDA.mat')
load('fEPSP_LFS1Hz_450p_SKF3m_m5.mat')
load('fEPSP_LFS1Hz_450p_SKF3m_m5_keinDA.mat')

t_1 = FT_LFS1Hz_450p_SKF3m_m5;
fEPSP_1 = fEPSP_LFS1Hz_450p_SKF3m_m5;
t_1_keinDA = FT_LFS1Hz_450p_SKF3m_m5_keinDA;
fEPSP_1_keinDA = fEPSP_LFS1Hz_450p_SKF3m_m5_keinDA;

load('chen_LFS.mat')
load('chen_LFS_SKF.mat')

t1 = chen_LFS(:,1);
x1 = chen_LFS(:,2);

t2 = chen_LFS_SKF(:,1);
x2 = chen_LFS_SKF(:,2);



ls = 1.5;
ms = 6;


a = 1;
b = 17;
c = 22;

HFS = 33*[1,1]-10; 
SKF = [0,15]+33-5-10;

figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
plot(t1+33-10,x1,'go',t2+23-10,x2,'ro',...
    t_1_keinDA(a:b)./60000-10,fEPSP_1_keinDA(a:b),'ks-', t_1(a:b)./60000-10,fEPSP_1(a:b),'b^-',...
    HFS,[105,115],'k-',SKF,[120,120],'b',...
    t_1_keinDA(c:end)./60000-10,fEPSP_1_keinDA(c:end),'ks-',t_1(c:end)./60000-10,fEPSP_1(c:end),'b^-',...
    'LineWidth',ls,'MarkerSize',ms)
ylabel('fEPSP (% change)')
xlabel('Time (min)')
xlim([0,80])
ylim([30,125])
legend('(Chen, 1995) - LFS','(Chen, 1995) - LFS + SKF','(Model) - LFS','(Model) - LFS + SKF')





