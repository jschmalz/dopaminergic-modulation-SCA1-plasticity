clearvars
clc
clear -all


load('FT_LFS1Hz_900b_SCH.mat')
load('FT_LFS1Hz_900b_SCH_keinDA.mat')
load('fEPSP_LFS1Hz_900b_SCH.mat')
load('fEPSP_LFS1Hz_900b_SCH_keinDA.mat')

t_1 = FT_LFS1Hz_900b_SCH;
fEPSP_1 = fEPSP_LFS1Hz_900b_SCH;
t_1_keinDA = FT_LFS1Hz_900b_SCH_keinDA;
fEPSP_1_keinDA = fEPSP_LFS1Hz_900b_SCH_keinDA;

load('frey_2004_SCH.mat')
load('Frey_2004_SLFS.mat')

xt = frey_2004_SCH(:,1);
xm = frey_2004_SCH(:,2);


ls = 1.5;
ms = 6;

HFS = 33*[1,1]; 
SKF = [0,20]+33+900*1000/60000;

a = 1;
b = 17;
c = 26;
npts = 5;

figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
plot(xt+31,xm,'ro',...
    t_1_keinDA(a:npts:b)./60000,fEPSP_1_keinDA(a:npts:b),'ks-', t_1(a:npts:b)./60000,fEPSP_1(a:npts:b),'b^-',...
    HFS,[105,115],'k-',...
    t_1_keinDA(c:npts:end)./60000,fEPSP_1_keinDA(c:npts:end),'ks-',t_1(c:npts:end)./60000,fEPSP_1(c:npts:end),'b^-',...
    'LineWidth',ls,'MarkerSize',ms)
ylabel('fEPSP (% change)')
xlabel('Time (min)')
xlim([0,550])
ylim([30,120])
legend('(Frey, 2004) - LFS + SCH','(Model) - LFS','(Model) - LFS + SCH')
xticks([0:100:600])
