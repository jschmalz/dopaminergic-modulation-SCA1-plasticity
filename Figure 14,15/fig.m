clearvars
clc
clear -all


load('FT_LFS3Hz_2400p_SKF100m_0.mat')
load('FT_LFS3Hz_2400p_SKF100m_0_keinDA.mat')
load('fEPSP_LFS3Hz_2400p_SKF100m_0.mat')
load('fEPSP_LFS3Hz_2400p_SKF100m_0_keinDA.mat')

load('FT_LFS3Hz_2x1200p_SKF100m_0.mat')
load('FT_LFS3Hz_2x1200p_SKF100m_0_keinDA.mat')
load('fEPSP_LFS3Hz_2x1200p_SKF100m_0_keinDA.mat')
load('fEPSP_LFS3Hz_2x1200p_SKF100m_0.mat')



t_1 = FT_LFS3Hz_2400p_SKF100m_0;
fEPSP_1 = fEPSP_LFS3Hz_2400p_SKF100m_0;
t_1_keinDA = FT_LFS3Hz_2400p_SKF100m_0_keinDA;
fEPSP_1_keinDA = fEPSP_LFS3Hz_2400p_SKF100m_0_keinDA;

t_2 = FT_LFS3Hz_2x1200p_SKF100m_0;
fEPSP_2 = fEPSP_LFS3Hz_2x1200p_SKF100m_0;
t_2_keinDA = FT_LFS3Hz_2x1200p_SKF100m_0_keinDA;
fEPSP_2_keinDA = fEPSP_LFS3Hz_2x1200p_SKF100m_0_keinDA;



load('mockett_2x_1200_SKF.mat')
load('mockett2400_SKF.mat')

td_m2x = mockett_2x_1200_SKF(:,1);
xd_m2x = mockett_2x_1200_SKF(:,2);

td_m = mockett2400_SKF(:,1);
xd_m = mockett2400_SKF(:,2);


load('fEPSP_LFS1Hz_900p_SKF100m_0.mat')
load('FT_LFS1Hz_900p_SKF100m_0.mat')
load('fEPSP_LFS1Hz_900p_SKF100m_0_keinDA.mat')
load('FT_LFS1Hz_900p_SKF100m_0_keinDA.mat')

fEPSP_3 = fEPSP_LFS1Hz_900p_SKF100m_0;
t_3 = FT_LFS1Hz_900p_SKF100m_0;
fEPSP_3_keinDA = fEPSP_LFS1Hz_900p_SKF100m_0_keinDA;
t_3_keinDA = FT_LFS1Hz_900p_SKF100m_0_keinDA;


load('fEPSP_LFS1Hz_900b_SKF100m_0.mat')
load('fEPSP_LFS1Hz_900b_SKF100m_0_keinDA.mat')
load('FT_LFS1Hz_900b_SKF100m_0.mat')
load('FT_LFS1Hz_900b_SKF100m_0_keinDA.mat')

fEPSP_4 = fEPSP_LFS1Hz_900b_SKF100m_0;
t_4 = FT_LFS1Hz_900b_SKF100m_0;
fEPSP_4_keinDA = fEPSP_LFS1Hz_900b_SKF100m_0_keinDA;
t_4_keinDA = FT_LFS1Hz_900b_SKF100m_0_keinDA;




ls = 1.5;
ms = 6;

HFS = 33*[1,1]; 
SKF = [0,20]+33+2*6;

nopts = 5;

N1 = length(t_1);
I1 = [1:nopts:N1];

a = 1;
b = 17;
c = 24;

 load('mocket_2xLFS.mat')
 load('mockett_2x_1200.mat')
 
 %% Figure 14

figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
subplot(2,1,1)
plot(mocket_2xLFS(:,1)+1,mocket_2xLFS(:,2)+100,'go',td_m+1,xd_m+100,'ro',...
    t_1_keinDA(a:b)./60000,fEPSP_1_keinDA(a:b),'ks-', t_1(a:b)./60000,fEPSP_1(a:b),'b^-',...
    HFS,[110,120],'k-',SKF,[115,115],'b',...
    t_1_keinDA(c:end)./60000,fEPSP_1_keinDA(c:end),'ks-',t_1(c:end)./60000,fEPSP_1(c:end),'b^-',...
    'LineWidth',ls,'MarkerSize',ms)
ylabel('fEPSP (% change)')
xlabel('Time (min)')
xlim([0,160])
ylim([40,130])
legend('(Mockett, 2007) - LFS','(Mockett, 2007) - LFS + SKF','(Model) - LFS','(Model) - LFS + SKF')




 %% Figure 15
 
a = 1;
b = 17;
c = 21;
d = 22;
e = 27;

HFS = 33*[1,1]; 
SKF = [0,20]+33+2*6+5;

subplot(2,1,2)
plot(mockett_2x_1200(:,1)+1,mockett_2x_1200(:,2)+100,'go',td_m2x+1,xd_m2x+100,'ro',...
     t_2_keinDA(a:b)./60000,fEPSP_2_keinDA(a:b),'ks-',t_2(a:b)./60000,fEPSP_2(a:b),'b^-',...
     HFS,[110,120],'k-',HFS+6+5,[110,120],'k-',SKF,[115,115],'b',...
     t_2_keinDA(c:d)./60000,fEPSP_2_keinDA(c:d),'ks-', t_2(c:d)./60000,fEPSP_2(c:d),'b^-',...
     t_2_keinDA(e:end)./60000,fEPSP_2_keinDA(e:end),'ks-',t_2(e:end)./60000,fEPSP_2(e:end),'b^-',...
     'LineWidth',ls,'MarkerSize',ms)
ylabel('fEPSP (% change)')
xlabel('Time (min)')
xlim([0,160])
ylim([40,130])
legend('(Mockett, 2007) - LFS','(Mockett, 2007) - LFS + SKF','(Model) - LFS','(Model) - LFS + SKF')



%%



a = 1;
b = 17;
c = 26;
d = 30;
e = 38;
f = 42;
g = 51;


HFS = 33*[1,1]; 
SKF = [0,20] + 33 + 3*900*1000/60000 + 2*10;
stim_duration = 15;
IBI = 10;

figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
subplot(2,1,1)
plot(t_3_keinDA(a:b)./60000,fEPSP_3_keinDA(a:b),'ks-',t_3(a:b)./60000,fEPSP_3(a:b),'b^-',...
     HFS,[110,120],'k-',HFS+IBI+stim_duration,[110,120],'k-',HFS+2*IBI+2*stim_duration,[110,120],'k-',SKF,[115,115],'b',...
     t_3_keinDA(c:d)./60000,fEPSP_3_keinDA(c:d),'ks-', t_3(c:d)./60000,fEPSP_3(c:d),'b^-',...
     t_3_keinDA(e:f)./60000,fEPSP_3_keinDA(e:f),'ks-',t_3(e:f)./60000,fEPSP_3(e:f),'b^-',...
     t_3_keinDA(g:end)./60000,fEPSP_3_keinDA(g:end),'ks-',t_3(g:end)./60000,fEPSP_3(g:end),'b^-',...
     'LineWidth',ls,'MarkerSize',ms)
ylabel('fEPSP (% change)')
xlabel('Time (min)')
ylim([30,125])
xlim([0,250])
legend('(Model) - LFS','(Model) - LFS + SKF')


a = 1;
b = 17;
c = 26;
d = 40;


HFS = 33*[1,1]; 
SKF = [0,20] + 33 + 1*900*1000/60000 + 0*10;

subplot(2,1,2)
plot(t_4_keinDA(a:b)./60000,fEPSP_4_keinDA(a:b),'ks-',t_4(a:b)./60000,fEPSP_4(a:b),'b^-',...
     HFS,[110,120],'k-',SKF,[115,115],'b',...
     t_4_keinDA(c:end)./60000,fEPSP_4_keinDA(c:end),'ks-', t_4(c:end)./60000,fEPSP_4(c:end),'b^-',...
     'LineWidth',ls,'MarkerSize',ms)
ylabel('fEPSP (% change)')
xlabel('Time (min)')
ylim([30,125])
xlim([0,250])
legend('(Model) - LFS','(Model) - LFS + SKF')


