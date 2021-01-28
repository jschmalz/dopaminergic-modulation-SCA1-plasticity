clearvars
clc
clear -all



load('zhang_2009_LTP.mat')
td_1 = zhang_2009_LTP(:,1)-37;
xd_1 = zhang_2009_LTP(:,2)*100;

load('blitzer_1995_LTP.mat')
td_4 = blitzer_1995_LTP(:,1)+5;
xd_4 = blitzer_1995_LTP(:,2);

% new blizter, 1995
bliz = [1.276769899629656, 211.64497343137776
5.059443937523483, 196.99801406258385
10.121571574257956, 185.1256642155547
15.388331275830609, 182.66404111427187
20.29346250872203, 181.03067468198162
29.93196822500134, 181.91428264720088
40.12532875315334, 186.6755005099028
50.13150126133864, 189.22166013633188
60.11955879984971, 184.29539477215388
69.93988513767378, 185.18000912457728
79.95612151790026, 191.87751596801027
89.75833288605013, 185.28970532982657
99.74907412377222, 181.470465890183];

td_4b = bliz(:,1);
xd_4b = bliz(:,2);


load('Stramiello_HFS.mat')
td_6 = Stramiello_HFS(:,1)-44+1;
xd_6 = Stramiello_HFS(:,2)*100;

load('Li_2013.mat')
td_8 = Li_2013(:,1);
xd_8 = Li_2013(:,2);

load('hernandez_HFS_300pulses.mat')
t_17 = hernandez_HFS_300pulses(:,1)+1.5;
x_17 = hernandez_HFS_300pulses(:,2);


load('Karpova_2006.mat')
td_9 = Karpova_2006(:,1)+2;
xd_9 = Karpova_2006(:,2);

load('kasashara_2001_LTP.mat')

td_10 = kasashara_2001_LTP(:,1)-20+1;
xd_10 = kasashara_2001_LTP(:,2);

load('roberto_2003_LTP.mat')

td_5 = roberto_2003_LTP(:,1);
xd_5 = roberto_2003_LTP(:,2);

load('hernandez_HFS_200pulses.mat')
t_16 = hernandez_HFS_200pulses(:,1)+1;
x_16 = hernandez_HFS_200pulses(:,2);

load('klar_2015_LTP.mat')
td_13 = klar_2015_LTP(:,1)+1;
xd_13 = klar_2015_LTP(:,2);

load('breakwell_1998_LTP.mat')
td_11 = breakwell_1998_LTP(:,1);
xd_11 = breakwell_1998_LTP(:,2)*100;

load('Papatheodoropoulos_1998.mat')
td_12 = Papatheodoropoulos_1998(:,1)-10;
xd_12 = Papatheodoropoulos_1998(:,2);

load('hernandez_HFS_100pulses.mat')
t_15 = hernandez_HFS_100pulses(:,1);
x_15 = hernandez_HFS_100pulses(:,2);

td_5 = roberto_2003_LTP(:,1)+2.5;
xd_5 = roberto_2003_LTP(:,2);

load('hernandez_HFS_200pulses.mat')
t_16 = hernandez_HFS_200pulses(:,1)+1;
x_16 = hernandez_HFS_200pulses(:,2);

load('klar_2015_LTP.mat')
td_13 = klar_2015_LTP(:,1)+1;
xd_13 = klar_2015_LTP(:,2);

load('breakwell_1998_LTP.mat')
td_11 = breakwell_1998_LTP(:,1);
xd_11 = breakwell_1998_LTP(:,2)*100;

load('Papatheodoropoulos_1998.mat')
td_12 = Papatheodoropoulos_1998(:,1)-10+1.1;
xd_12 = Papatheodoropoulos_1998(:,2);

load('hernandez_HFS_100pulses.mat')
t_15 = hernandez_HFS_100pulses(:,1)+2;
x_15 = hernandez_HFS_100pulses(:,2);

%%
load('fEPSP_3xHFS_10min.mat')
load('t_3xHFS_10min.mat')

load('fEPSP_3xHFS_0p5s.mat')
load('t_3xHFS_0p5s.mat')

load('fEPSP_3xHFS_20s.mat')
load('t_3xHFS_20s.mat')

load('fEPSP_3xHFS_10s.mat')
load('t_3xHFS_10s.mat')


ls = 1.5;
ms = 6;

% figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
subplot(2,2,1)
plot(td_1,xd_1,'ro',t_3xHFS_0p5s/60000-10,fEPSP_3xHFS_0p5s,'ks-','LineWidth',ls,'MarkerSize',ms)
legend('(Zhang, 2009)', '(Model)')
xlabel('Time (min)')
ylabel('fEPSP (% change)')
xlim([-5 50])
ylim([95 300])
title('3x 100 Hz for 1s w/ 0.5 s interval')

s = 20 + 5/60;
subplot(2,2,2)
plot(td_4b+s,xd_4b,'ro',td_9-1,xd_9,'bo',t_3xHFS_10min/60000-10,fEPSP_3xHFS_10min,'ks-','LineWidth',ls,'MarkerSize',ms)
legend('(Blitzer, 1995)','(Karpova, 2006)','(Model)')
xlabel('Time (min)')
ylabel('fEPSP (% change)')
xlim([-5,110])
ylim([95 300])
title('3x 100 Hz for 1s w/ 10 min interval')

subplot(2,2,3)
plot(td_6-1,xd_6,'ro',t_3xHFS_20s/60000-10,fEPSP_3xHFS_20s,'ks-','LineWidth',ls,'MarkerSize',ms)
legend('(Stramuiello, 2008)', '(Model)')
xlabel('Time (min)')
ylabel('fEPSP (% change)')
xlim([-5,30])
ylim([95 300])
title('3x 100 Hz for 1 sec w/ 20 sec interval')

subplot(2,2,4)
plot(t_17-1,x_17,'ro',t_3xHFS_10s/60000-10,fEPSP_3xHFS_10s,'ks-','LineWidth',ls,'MarkerSize',ms)
legend('(Hernandez, 2005)','(Model)')
xlabel('Time (min)')
ylabel('fEPSP (% change)')
xlim([-5,60])
ylim([95 300])
title('3x 100 Hz for 1 sec w/ 10 sec interval')

%%

load('fEPSP_2xHFS_10s.mat')
load('t_2xHFS_10s.mat')

load('fEPSP_2xHFS_20s.mat')
load('t_2xHFS_20s.mat')


figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')

subplot(2,2,2)
plot(td_10-1,xd_10,'ro',t_2xHFS_20s/60000-10,fEPSP_2xHFS_20s,'ks-','LineWidth',ls,'MarkerSize',ms)
legend('(Kasashara, 2001)')
xlabel('Time (min)')
ylabel('fEPSP (% change)')
xlim([-5,60])
ylim([95 220])
title('2x 100 Hz for 1 sec w/ 20 sec interval')

subplot(2,2,3)
plot(t_16-1,x_16,'ro',t_2xHFS_10s/60000-10,fEPSP_2xHFS_10s,'ks-','LineWidth',ls,'MarkerSize',ms)
legend('(Hernandez, 2005)')
xlabel('Time (min)')
ylabel('fEPSP (% change)')
xlim([-5,60])
ylim([95 220])
title('2x 100 Hz for 1 sec w/ 10 sec interval')



load('fEPSP_4xHFS_5min.mat')
load('t_4xHFS_5min.mat')

load('fEPSP_1xHFS.mat')
load('t_1xHFS.mat')

% figure(3)
subplot(2,2,1)
plot(td_5-1,xd_5,'ro',td_12-1,xd_12,'bo',t_15-1,x_15,'go',t_1xHFS/60000-10,fEPSP_1xHFS,'ks-','LineWidth',ls,'MarkerSize',ms)
legend('(Roberto, 2003)','(Papatheodoropoulos, 1998)','(Hernandez, 2005)')
xlabel('Time (min)')
ylabel('fEPSP (% change)')
ylim([95 220])
xlim([-5,60])
title('1 HFS train')

subplot(2,2,4)
plot(td_8,xd_8,'ro',t_4xHFS_5min/60000-10,fEPSP_4xHFS_5min,'ks-','LineWidth',ls,'MarkerSize',ms)
legend('(Li, 2013)','(Model)')
xlabel('Time (min)')
ylabel('fEPSP (% change)')
ylim([95 300])
xlim([-5,110])
title('4x HFS trains with 5 minute interval')

