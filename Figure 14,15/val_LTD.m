clear -all
clc
clearvars

load('theta_LTD_08192020.mat')
theta_P = mean(theta_LTD_08192020,2);

gamma = theta_P(1);
eta = theta_P(2);
nu = theta_P(3);
g = theta_P(4);
M_d = theta_P(5);
p_d = theta_P(6);
A_d = theta_P(7);
M_p = theta_P(8);
p_p = theta_P(9);
A_p = theta_P(10);
f = 281.3122;

% equilibrium parameters
p_plast(1) = f;
p_plast(2) = M_d;
p_plast(3) = p_d;
p_plast(4) = A_d;
p_plast(5) = M_p;
p_plast(6) = p_p;
p_plast(7) = A_p;


gamma = theta_P(1);
eta = theta_P(2);
nu = theta_P(3);
g = theta_P(4);


% drug properties
m = 1;
EC50_SKF = 20e-6;
EC50_APB = 15e-6;
EC50_DA = 40e-06;

Emax_SKF = 0.6;
Emax_APB = 0.4;
Emax_DA = 0.6;

nH_SKF = 1;
nH_APB = 1;
nH_DA = 2;

IC50_SKF = 1e-07;
IC50_DA = 0.1e-06;
IC50_APB = 0.1e-07;

Imax_SKF = -0.05; 
Imax_APB = -0.03;
Imax_DA = -0.25;

DP(1) = EC50_SKF;
DP(2) = IC50_SKF;
DP(3) = Emax_SKF;
DP(4) = Imax_SKF;
DP(5) = nH_SKF  ;



duration = 200*60*1000;
dstim = 300;
freq = 2*60*1000;
start_stim = 33*60*1000;
stim_duration = 1200*300;
three = 0;
IBI = 5*60*1000;
four = 0;
two = 0;
start_SKF = 33*60*1000 + 1*6.7*60*1000 + 0*60*60*1000 + 0*IBI;
three_SKF = 0;
SKF_duration = 20;
LTD = 1;
LFS_prot = 0;
SKF_conc = 100e-6;
basal_DA = 1;

P_run(1) = duration;
P_run(2) = dstim;
P_run(3) = freq;
P_run(4) = start_stim;
P_run(5) = stim_duration;
P_run(6) = three;
P_run(7) = IBI;
P_run(8) = four;
P_run(9) = two;
P_run(10) = start_SKF;
P_run(11) = three_SKF;
P_run(12) = SKF_duration;
P_run(13) = SKF_conc;
P_run(14) = LTD;
P_run(15) = LFS_prot;
P_run(16) = basal_DA;






disp('0/2')
tic
%run model EPSP predict
[fEPSP, Ft, LTP, time,TLTP_ON,TLTD_ON,v,tag_DA,tag_stim] = CA1_model(theta_P,P_run,p_plast,DP);
disp('Run Time (min):')
disp(toc/60)


figure(3)
subplot(2,1,1)
plot(time/6000,tag_DA,'k')

subplot(2,1,2)
plot(time/6000,tag_stim,'b')


SKF_conc = 0;
basal_DA = 1;
P_run(13) = SKF_conc;
P_run(16) = basal_DA;



disp('1/2')
tic
%run model EPSP predict
[fEPSP_keinDA, Ft_keinDA, LTP_kda, time_kda,TLTP_ON_kda,TLTD_ON_kda,v_kda,tag_DAkda,tag_stimkda] = CA1_model(theta_P,P_run,p_plast,DP);
disp('Run Time (min):')
disp(toc/60)





figure(1)
subplot(2,1,1)
plot(Ft_keinDA/(60*1000),fEPSP_keinDA,'ks',Ft/(60*1000),fEPSP,'bs')
xlabel('Time (min)')
ylabel('fEPSP (% change)')

subplot(2,1,2)
plot(Ft_keinDA./(60*1000),fEPSP-fEPSP_keinDA,'ks-')
xlabel('Time (min)')
ylabel('\Delta fEPSP')



% 
% % 
% % 
fEPSP_LFS3Hz_1200p_SKF100m_m0 = fEPSP;
FT_LFS3Hz_1200p_SKF100m_m0 = Ft;
save('fEPSP_LFS3Hz_1200p_SKF100m_m0.mat','fEPSP_LFS3Hz_1200p_SKF100m_m0')
save('FT_LFS3Hz_1200p_SKF100m_m0.mat','FT_LFS3Hz_1200p_SKF100m_m0')



fEPSP_LFS3Hz_1200p_SKF100m_m0_keinDA = fEPSP_keinDA;
FT_LFS3Hz_1200p_SKF100m_m0_keinDA = Ft_keinDA;
save('fEPSP_LFS3Hz_1200p_SKF100m_m0_keinDA.mat','fEPSP_LFS3Hz_1200p_SKF100m_m0_keinDA')
save('FT_LFS3Hz_1200p_SKF100m_m0_keinDA.mat','FT_LFS3Hz_1200p_SKF100m_m0_keinDA')

% [slope_EPSP2, F_t2] = make_slopeEPSP(v_kda,time_kda,freq)
