%HH ODE
function dxdt = HHODE_AMPA_NMDA(t,x,xP,xP2,xP3,xP4,stim_ON)
% state
V = x(1,1);
hNaT = x(2,1);
nKDR = x(3,1);
bKA = x(4,1);
uKM = x(5,1);
rCa = x(6,1);
Cai = x(7,1);
cCa = x(8,1);
qCa = x(9,1);
h_AMPA = x(10,1);
U_AMPA = x(11,1);
h_NMDA = x(12,1);
U_NMDA = x(13,1);
C = x(14,1);
Np = x(15,1);
Nd = x(16,1);

TPKA_ON = x(17,1);
TPKA_OFF = x(18,1);
TPLC_ON = x(19,1);
TPLC_OFF = x(20,1);
TLTP_ON = x(21,1);
TLTD_ON = x(22,1);
tag_stim_off = x(23,1);
tag_stim = x(24,1);
tag_DA_off = x(25,1);
tag_DA = x(26,1);
T_OFF = x(27,1);
HFS = x(28,1);



% neuron membrane parameters
Cm = xP(1,1);
I = xP(2,1);
gNaT = xP(3,1);
ENa = xP(4,1);
gNaP = xP(5,1);
gKDR = xP(6,1);
EK = xP(7,1);
gKA = xP(8,1);
gKM = xP(9,1);
EL = xP(10,1);
gL = xP(11,1);
gCa = xP(12,1);
gC =  xP(13,1); 
gsAHP = xP(14,1); 
ECa = xP(15,1);


% AMPA and NMDA paramters

alpha = xP2(1);
Es = xP2(2); %mV
tau1_ampa = xP2(3); %ms
tau2_ampa = xP2(4); %ms
gsbar = xP2(5); %nS
tau1_nmda = xP2(6); %ms
tau2_nmda = xP2(7); %ms
mu = xP2(8); %mM-1
Mg = xP2(9); % mM
gamma_nmda = xP2(10); %mV-1


% LTP - LTD parameters
tot_Ci = xP3(14);
rhoC = 1 + 0.31*(tot_Ci>3);
rhoB =  1 - 0.65*(tot_Ci>5);

gamma = xP3(1); %ms-1
eta = xP3(2);
nu = xP3(3);
g = xP3(4); %nA-1ms-1
f = xP3(5); %nA-1
p_d = rhoB*xP3(6); %ms-1
p_p = rhoC*xP3(7); %ms-1
M_p = rhoC*xP3(8); %nA/msec
M_d = rhoB*xP3(9); %nA/ms
A_p = rhoC*xP3(10); %nA2
A_d = rhoB*xP3(11); %nA2
C_m_p = xP3(12); % pF
n_syn = xP3(13);



% DA dynamic parameters
k1 =   8.2575e-06;
k2 =   7.1029e-08;
k3 =   3.2390e-06;
k4 =   1.4430e-07;
k5 =   2.8866e-07;
k6 =   2.3514e-06;
k7 =   1.6253e-07;
k8 =   5.3205e-04;

hillE = xP4(1);
hillI = xP4(2);
Basal_DA = xP4(3);
E0 = xP4(4);
I0 = xP4(5);
L_ON = xP4(6);

% stim + DAergic
kHFS =   1.3444e-07;
kmodstim =  556.7311;
knonlin =   1.5143e-06;
kLTD =   4.3622e-07;
knonlin_I =    2.4902;
tau_e =   2.4563e+05;
tau_d =   7.7315e+04;
kDA =    0.0010;
kSOP = kHFS;


fun_gstofEPSP = @(x) 0.625414*(x) + 37.3962;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plasticity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gmax_NMDA = (1).*gsbar./(1 + mu.*Mg.*exp(-gamma_nmda.*V));
g_NMDA = gmax_NMDA*h_NMDA;

I_caNMDA = 0.1*alpha*g_NMDA*(V-20);

dCdt =  (-gamma*stim_ON*I_caNMDA - eta*C); 
dNpdt = (nu*C -  (p_p - stim_ON*I_caNMDA*g)*Np +  (M_p*Np^2)/(A_p+Np^2));
dNddt = (nu*C - (p_d - stim_ON*I_caNMDA*g)*Nd + (M_d*Nd^2)/(A_d+Nd^2)); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DAerig and DA + Stim interactions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dTPKA_ONdt =  k1*hillE*TPKA_OFF - k2*TPKA_ON - k7*TPKA_ON*TPLC_ON;
dTPKA_OFFdt = -dTPKA_ONdt;
dTPLC_ONdt =  k3*abs(hillI)*TPLC_OFF - k4*TPLC_ON - k8*TPLC_ON*hillE^2 ; %
dTPLC_OFFdt = -dTPLC_ONdt;
dTLTP_ONdt = k5*TPKA_ON.^2*T_OFF;
dTLTD_ONdt = k6*TPLC_ON.^2*T_OFF;


dtag_stim_offdt = -kmodstim*abs(dNpdt-dNddt)*(abs(dNpdt-dNddt)>2e-8)*tag_stim_off + tag_stim./tau_d;
dtag_stimdt =  kmodstim*abs(dNpdt-dNddt)*(abs(dNpdt-dNddt)>2e-8)*tag_stim_off - tag_stim./tau_d;
dtag_DA_offdt = -kDA*abs(hillE-knonlin_I*abs(hillI))*tag_DA_off + tag_DA./tau_e;
dtag_DAdt = kDA*abs(hillE-knonlin_I*abs(hillI))*tag_DA_off - tag_DA./tau_e;


%%%%%

dT_OFFdt = -dTLTD_ONdt  - dTLTP_ONdt - kLTD*T_OFF.*(abs(Np-Nd)>10^-4)...
    - sign(E0-knonlin_I*abs(I0))*knonlin*abs(tag_stim*tag_DA) ; 

dHFSdt = f*(dNpdt-dNddt) - kHFS*HFS*(Basal_DA == 0) ... 
    + (dTLTP_ONdt - dTLTD_ONdt)...
    + sign(E0-knonlin_I*abs(I0))*knonlin*tag_DA*tag_stim ...
    - kSOP*(-1.12 + HFS)*(L_ON)*(fun_gstofEPSP(HFS*100+100) > 170) ;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AMPA and NMDA currents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dh_AMPAdt = U_AMPA;
dU_AMPAdt = -U_AMPA/tau2_ampa - (U_AMPA/tau1_ampa) - (h_AMPA/(tau1_ampa*tau2_ampa));
dh_NMDAdt = U_NMDA;
dU_NMDAdt = -U_NMDA/tau2_nmda - (U_NMDA/tau1_nmda) - (h_NMDA/(tau1_nmda*tau2_nmda));

gmax_AMPA = (1 + HFS).*gsbar.*tau1_ampa.*tau2_ampa./(tau2_ampa-tau1_ampa);
g_AMPA = gmax_AMPA*h_AMPA;


gmax_NMDA = (1).*gsbar./(1 + mu.*Mg.*exp(-gamma_nmda.*V));
g_NMDA = gmax_NMDA*h_NMDA;

gs = (1-alpha)*g_AMPA + alpha*g_NMDA;
Isyn =  n_syn*gs*(V-Es); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Na currents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Na Transient Current
thetam = -30;
sigm = 9.5;

thetah = -45;
sigh = -7;

thetaht = -40.5;
sight = -6;
mNaTinf = 1/(1+exp(-(V-thetam)/sigm));
tauhNaT = 1.0+7.5/(1+exp(-(V-thetaht)/sight));

hNaTinf = 1/(1+exp(-(V-thetah)/sigh));
dhNaT = 10*(hNaTinf-hNaT)/tauhNaT;
INaT = gNaT*mNaTinf^3*hNaT*(V-ENa);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Persistent Na Current
thetap = -46; 
sigp = 3;
mNaPinf = 1/(1+exp(-(V-thetap)/sigp));
 
INaP = gNaP*mNaPinf*(V-ENa);



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %K currents
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Delayed rectifier K current
thetan = -35;
signo = 10;
thetant = -27;
signt = -15;

nKDRinf = 1/(1+exp(-(V-thetan)/signo));
taunKDR = 1+5/(1+exp(-(V-thetant)/signt));
dnKDR = 10*(nKDRinf-nKDR)/taunKDR;
IKDR = gKDR*nKDR^4*(V-EK);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A-type K current
thetaa = -50;
siga = 20;
thetab = -80;
sigb = -6;

aKAinf = 1/(1+exp(-(V-thetaa)/siga));



%b inactivation
bKAinf = 1/(1+exp(-(V-thetab)/sigb));
taubKA = 15;

dbKA = (bKAinf-bKA)/taubKA;
IKA = gKA*aKAinf^3*bKA*(V-EK);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%M-type K current
thetaz = -39;
sigz = 5;

uKMinf = 1/(1+exp(-(V-thetaz)/sigz));
tauuKM =75;

duKM = (uKMinf-uKM)/tauuKM;
IKM = gKM*uKM*(V-EK);

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Leak current
IL = gL*(V-EL);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%High Threshold Ca current
thetar = -20;
sigr = 10;

taur = 1;
rinf = 1/(1+exp(-(V-thetar)/sigr));

drCa = (rinf-rCa)/taur;
ICa = gCa*rCa^2*(V-ECa);

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Ca dynamics
vCa = 0.13;
tauCa = 13;
dCaidt = -vCa*ICa - Cai/tauCa;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fast Ca activated K current
thetac = -30;
sigc = 7;

tauc = 2;
cinf = 1/(1+exp(-(V-thetac)/sigc));

ac = 6;
dinf = 1/(1+(ac/Cai));
dcCa = (cinf-cCa)/tauc;
IC = gC*dinf*cCa*(V-EK);

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Slow Ca activated K current

tauq = 450;
aq = 2;
qinf = 1/(1+(aq^4/Cai^4));
dqCa = (qinf-qCa)/tauq;
IsAHP = gsAHP*qCa*(V-EK);

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dVdt = (1/Cm)*(I -INaT - INaP - IKDR - IKA - IKM - IL - ICa - IsAHP -IC - Isyn ); %I is the external current



dxdt = [dVdt,dhNaT,dnKDR,dbKA,duKM,drCa,dCaidt,dcCa,dqCa,dh_AMPAdt,dU_AMPAdt,dh_NMDAdt,dU_NMDAdt,dCdt,dNpdt,dNddt,...
    dTPKA_ONdt,dTPKA_OFFdt,dTPLC_ONdt,dTPLC_OFFdt,dTLTP_ONdt,dTLTD_ONdt,dtag_stim_offdt,dtag_stimdt,dtag_DA_offdt...
    dtag_DAdt,dT_OFFdt,dHFSdt]';
end

