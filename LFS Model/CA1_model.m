function [fEPSP, F_t, LTP, time,TLTP_ON,TLTD_ON,v,tag_DA,tag_stim] = CA1_model(theta_P,P_run,p_plast,DP)
% run parameters
duration = P_run(1);
dstim = P_run(2);
freq = P_run(3);
start_stim = P_run(4);
stim_duration = P_run(5);
three = P_run(6);
IBI = P_run(7);
four = P_run(8);
two = P_run(9);
start_SKF = P_run(10);
three_SKF = P_run(11);
SKF_duration = P_run(12);
L = P_run(13);
LTD = P_run(14);
LFS_prot = P_run(15);
basal_DA = P_run(16);

first_one = 1;

%Model Paramters
C = 1; %mF/cm^2
I = 0;
gNaT = 35; %mS/cm^2
ENa = 55; %mV
gNaP = 0.3; %mS/cm^2 % Varies between 0 and 0.41 %0.2 %0.3
gKDR = 6; %mS/cm^2
EK = -90; %mV
gKA = 1.4; %mS/cm^2
gKM = 1; %mS/cm^2
EL = -70; %mV
gL = 0.05; %mS/cm^2
gCa = 0.08; %0.08 %0.02
gC = 10;
ECa = 120;
gsAHP = 5;
xP = [C, I, gNaT,ENa,gNaP,gKDR,EK,gKA,gKM,EL,gL,gCa,gC,gsAHP,ECa]';


% AMPA and NMDA paramters
alpha = 0.1;
Es = 0; %mV
tau1_ampa = 0.8263; %ms
tau2_ampa = 4.548; %ms
gsbar = 4.981e-6; %nS -> mS
tau1_nmda = 0.8189; %ms
tau2_nmda = 74.7877; %ms
mew_nmda = 0.2866; %mM-1
C_Mg = 1; % mM
gamma_nmda = 0.013; %mV-1
Uampa_stim = 0;
Unmda_stim = 0;
xP2 = [alpha,Es,tau1_ampa,tau2_ampa,gsbar,tau1_nmda,tau2_nmda,mew_nmda,C_Mg,gamma_nmda,Uampa_stim,Unmda_stim];


% LTP - LTD parameters
gamma = theta_P(1); %ms-1 
eta = theta_P(2);
nu = theta_P(3);
g = theta_P(4); %nA-1ms-1


f = p_plast(1); %nA-1
M_p = p_plast(2); %ms-1
p_p = p_plast(3); %ms-1
A_p = p_plast(4); %nA/msec 
M_d = p_plast(5); %nA/ms
p_d = p_plast(6); %nA2
A_d = p_plast(7); %nA2
C_m_p = 3; % pF
n_syn = 1.799e4;

totCi = 0;

xP3 = [gamma,eta,nu,g,f,p_d,p_p,M_p,M_d,A_p,A_d,C_m_p,n_syn,totCi];





% Initial Conditions
V0 = -71.81327;
hNaT0 = 0.98786;
nKDR0 = 0.02457;
bKA0 = 0.203517;
uKM0 = 0.00141;
rCa0 = 0.005507;
cCa0 = 0.002486;
qCa0 = 0;
Cai0 = 0.000787;
gampa0 = 0;
Uampa0 = 0;
Unmda0 = 0;
gnmda0 = 0;
C0 = 0;
Np0 = 0;
Nd0 = 0;
TPKA_ON0 = 0;
TPKA_OFF0 = 1;
TPLC_ON0 = 0;
TPLC_OFF0 = 1;
TLTP_ON0 = 0;
TLTD_ON0 = 0;
tag_stim_off0 = 1;
tag_stim0 = 0;
tag_DA_off0 = 1;
tag_DA0 = 0;
T_OFF0 = 1;
HFS0 = 0;

x0 = [V0,hNaT0,nKDR0,bKA0,uKM0,rCa0,Cai0,cCa0,qCa0,gampa0,Uampa0,gnmda0,Unmda0,C0,Np0,Nd0,...
    TPKA_ON0,TPKA_OFF0,TPLC_ON0,TPLC_OFF0,TLTP_ON0,TLTD_ON0,tag_stim_off0,tag_stim0,tag_DA_off0,tag_DA0,T_OFF0,HFS0]'; %Initial condition 





EC50 = DP(1);
IC50 = DP(2);
Emax = DP(3);
Imax = DP(4);
nH = DP(5);


fun_hill = @(EC50,Emax,L,nH) Emax./(1+(EC50./L).^nH);
E0 = fun_hill(EC50,Emax,L,nH);
I0 = - fun_hill(IC50,Imax,L,nH);


LFS_count=0;
LFS_ON = 0;
% run
for i = 1:ceil(duration/dstim)
    
    %ODE Call
    dt = 1;
    tspan = [dt:dt:dstim+dt];
    
    % --- parameter entry ---
    injection_time = start_SKF/(60*1000);
    hillE = E0*(i*dstim>injection_time*1000*60-1)*(i*dstim<(injection_time+SKF_duration)*1000*60) + E0*three_SKF*(i*dstim>(injection_time+2*SKF_duration)*1000*60-1)*(i*dstim<(injection_time+3*SKF_duration)*1000*60) + E0*three_SKF*(i*dstim>(injection_time+4*SKF_duration)*1000*60-1)*(i*dstim<(injection_time+5*SKF_duration)*1000*60);
    hillI = I0*(i*dstim>injection_time*1000*60-1)*(i*dstim<(injection_time+SKF_duration)*1000*60) + I0*three_SKF*(i*dstim>(injection_time+2*SKF_duration)*1000*60-1)*(i*dstim<(injection_time+3*SKF_duration)*1000*60) + I0*three_SKF*(i*dstim>(injection_time+4*SKF_duration)*1000*60-1)*(i*dstim<(injection_time+5*SKF_duration)*1000*60);

    xP4 = [hillE,hillI,basal_DA,E0,I0,(i*dstim>=start_SKF)]; 
    xP3 = [gamma,eta,nu,g,f,p_d,p_p,M_p,M_d,A_p,A_d,C_m_p,n_syn,totCi];
    
    % spike or EPSP generation
    EPSP_on = 0;
    if mod(dstim*i,freq) == 0 || dstim*i == 60*1000
       EPSP_on = 1; 
       first_one = 0;
    end
    


    if LFS_prot == 2
        AP_ON =(LFS_count>=0)*(LFS_count<150)*LFS_ON;         
    else
        AP_ON = (i*dstim>=start_stim)*(i*dstim<start_stim+stim_duration+dstim) + two*(i*dstim>=start_stim+stim_duration+IBI)*(i*dstim<start_stim+IBI+2*stim_duration+dstim) + three*(i*dstim>=start_stim+2*stim_duration+2*IBI)*(i*dstim<start_stim+2*IBI+3*stim_duration+dstim) + four*(i*dstim>=start_stim+3*IBI+3*stim_duration)*(i*dstim<start_stim+3*IBI+4*stim_duration+dstim);
    end
    EPSP = 0.42; 

    if (i*dstim >= start_stim)*(i*dstim <= start_stim+stim_duration)*(LFS_prot==2) == 1
        LFS_count = LFS_count + dstim;
        LFS_ON = 1;
        if LFS_count >= 1000-0.1
            LFS_count = 0;
        end
    else
        LFS_ON = 0;
    end
    
    stim_length = (IBI*(1+two+three+four) + stim_duration*(1+two+three+four));
    stim_ON = (i*dstim>=start_stim)*( i*dstim<=start_stim+stim_length);

    % ODE
    [t,x] = ode45(@(t,x) HHODE_AMPA_NMDA(t,x,xP,xP2,xP3,xP4,stim_ON), tspan, x0);
    
    % new ICs
    a = dstim*(i-1)/dt + 1;
    b = dstim*i/dt;
    c = size(x,1);
    uAMPA_0 = tau1_ampa^-1-tau2_ampa^-1;
    uNMDA_0 = tau1_nmda^-1-tau2_nmda^-1;
    x0 = x(c,:)';
    x0(11) = x(c,11) + EPSP_on*EPSP*uAMPA_0 + 1*AP_ON*uAMPA_0;
    x0(13) = x(c,13) + EPSP_on*EPSP*uNMDA_0 + 1*AP_ON*uNMDA_0;
    
    % record state variables
    v(a:b,1) = x(1:c-1,1);
    hampa(a:b,1) = x(1:c-1,10);
    hnmda(a:b,1) = x(1:c-1,12);
    Ci(a:b,1) = x(1:c-1,14);
    Np(a:b,1) = x(1:c-1,15);
    Nd(a:b,1) = x(1:c-1,16);
    Cai(a:b,1) = x(1:c-1,7);  
    HFS(a:b,1) = x(1:c-1,28); 
    TLTP_ON(a:b,1) = x(1:c-1,21);
    TLTD_ON(a:b,1) = x(1:c-1,22);
    tag_stim(a:b,1) = x(1:c-1,24); 
    tag_DA(a:b,1) = x(1:c-1,26); 
    
    if i>1
        time(a:b,1) = t(1:c-1,1) + time(a-1,1);
    else
        time(a:b,1) = t(1:c-1,1);
    end
   
    totCi = sum(Ci);
end

fun_gstofEPSP = @(x) 0.625414*(x) + 37.3962;


LTP = fun_gstofEPSP((HFS+1)*100);
LTP_raw = f.*(Np-Nd);
% convert to fEPSP
[slope_EPSP,F_t] = make_slopeEPSP(v,time,freq);
fEPSP_fun = @(x) 0.35912*x + 0.13199;
fEPSP = 100*fEPSP_fun(slope_EPSP)./fEPSP_fun(slope_EPSP(1));

totCi
end


