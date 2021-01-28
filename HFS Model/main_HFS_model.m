clear -all
clc
clearvars





%%%%%%%  SIMULATION PARAMETERS   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% change these parameters only
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
duration = 400*60*1000; % duration of simulation
dstim = 10;  % interval of HFS spikes 1/freq in ms
freq = 2*60*1000; % frequency of EPSP

% HFS parameters
start_stim = 214*60*1000 + 2000 ; % time of HFS delivery - little bit added to stimulation so it comes at different time then EPSP
stim_duration = 1000; % length of HFS in ms

% intertrain interval for HFS
IBI = 10*60*1000;  

% number of trains 2,3, or 4 (i.e for 4 trains two=1, three =1, and four =1)
two = 1;
three = 1;
four = 0;
ten = 0;

% DA agonist parameters
start_SKF = 2*60*1000 + 2000; % time of DA agonist injection
three_SKF = 0;  % if three DA agonist pulse three_SKF =1 else three_SKF = 0
SKF_duration = 15;  % duration of DA agonist pulse
SKF_conc = 50e-6; % DA agonist pulse

% SCH parameters
basal_DA = 1; % basal dopmine? if SCH present set basal_DA = 0


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


% changes these for agonist applied
DP(1) = EC50_SKF; 
DP(2) = IC50_SKF;
DP(3) = Emax_SKF;
DP(4) = Imax_SKF;
DP(5) = nH_SKF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
P_run(14) = ten;
P_run(15) = basal_DA;




% HFS parameters
gamma =    0.0364;
eta =    0.0035;
nu =    0.1070;
g =  214.1331;
M_p =   4.8347e-08;
p_p =   1.0066e-05;
A_p =   2.3018e-06;
M_d =   4.9091e-07;
p_d =   2.2314e-04;
A_d =  6.8059e-07;
f = 460.1983;



% % hand-fit to otm data for Figure 5 and 6
% gamma =  0.0364;
% eta =   0.0035;
% nu = 0.1070;
% g =  214.1331;
% M_p = 0.155*4.8347e-08;
% p_p = 0.15*1.0066e-05;
% A_p = 2.3018e-06;
% M_d = 4.9091e-07;
% p_d = 2.2314e-04;
% A_d = 6.8059e-07;
% f =  298;

p_plast(1) = f;
p_plast(2) = M_p;
p_plast(3) = p_p;
p_plast(4) = A_p;
p_plast(5) = M_d;
p_plast(6) = p_d;
p_plast(7) = A_d;

% run 
theta_P(1) = gamma;
theta_P(2) = eta;
theta_P(3) = nu;
theta_P(4) = g;







% % run parameters otm for figures 5 and 6
% duration = 150*60*1000;
% dstim = 10;
% freq = 60*500;
% start_stim = 10*60*1000 + 2000 ; % little bit added to stimulation comes at different time then EPSP
% stim_duration = 40;
% three = 1;
% IBI = 30;
% four = 1;  
% two = 1;
% start_SKF = 45*60*1000 + 2000;
% three_SKF = 0;
% SKF_duration = 5;
% SKF_conc = 5e-6;
% ten = 1;
% basal_DA = 1;







disp('0/2')
tic
%run model EPSP predict
[fEPSP,Ft,v,time,Ci,hampa,hnmda,Np,Nd] = CA1_model(theta_P,P_run,p_plast,DP);
disp('Run Time (min): ')
disp(toc/60)



disp('1/2')
SKF_conc = 0;
P_run(13) = SKF_conc;
P_run(15) = 1;
tic
[fEPSP_keinDA, Ft_keinDA] = CA1_model(theta_P,P_run,p_plast,DP);
disp('Run Time (min): ')
disp(toc/60)





figure(1)
subplot(2,1,1)
plot(Ft_keinDA/(60*1000),fEPSP_keinDA,'ks-',Ft/(60*1000),fEPSP,'bs-')
xlabel('Time (min)')
ylabel('fEPSP (% change)')

subplot(2,1,2)
plot(Ft_keinDA./(60*1000),fEPSP-fEPSP_keinDA,'ks-')
xlabel('Time (min)')
ylabel('\Delta fEPSP')



