clearvars
clc
clear -all

%% load parameter distributions

load('theta_SOP_08182020.mat')
load('theta_nov_24_2020.mat')
load('theta_LTD_08192020.mat')
load('theta_HFS_08192020.mat')


%% Figure 22
clc

theta_HFS = mean(theta_HFS_08192020,2);

f = 460.1983;
gamma = theta_HFS(1);
eta = theta_HFS(2);
nu = theta_HFS(3);
g = theta_HFS(4);
M_p = theta_HFS(5);
p_p = theta_HFS(6);
A_p = theta_HFS(7);
M_d = theta_HFS(8);
p_d = theta_HFS(9);
A_d = theta_HFS(10);


std_gamma = std(theta_HFS_08192020(1,:));
std_eta = std(theta_HFS_08192020(2,:));
std_nu = std(theta_HFS_08192020(3,:));
std_g = std(theta_HFS_08192020(4,:));
std_M_p = std(theta_HFS_08192020(5,:));
std_p_p = std(theta_HFS_08192020(6,:));
std_A_p = std(theta_HFS_08192020(7,:));
std_M_d = std(theta_HFS_08192020(8,:));
std_p_d = std(theta_HFS_08192020(9,:));
std_A_d = std(theta_HFS_08192020(10,:));

% %%
figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
subplot(2,5,1)
hold on;
line([gamma,gamma], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_HFS_08192020(1,:),'FaceColor','none','LineWidth',1.2)
ylabel('Number of Observations')
xlabel('\gamma_C (ms^{-1})')
hold on;
line([gamma,gamma], [0, 160], 'LineWidth', 2, 'Color', 'r');
legend(string(gamma)+' ms^{-1}')
ylim([0,160])

subplot(2,5,2)
hold on;
line([eta,eta], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_HFS_08192020(2,:),'FaceColor','none','LineWidth',1.2)
xlabel('\eta (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([eta,eta], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend(string(eta)+' ms^{-1}')

subplot(2,5,6)
hold on;
line([nu,nu], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_HFS_08192020(3,:),'FaceColor','none','LineWidth',1.2)
xlabel('\nu (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([nu,nu], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend(string(nu)+' ms^{-1}')

subplot(2,5,7)
hold on;
line([g,g], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_HFS_08192020(4,:),'FaceColor','none','LineWidth',1.2)
xlabel('g (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([g,g], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend(string(g)+' ms^{-1}')

% figure(2)
subplot(2,5,3)
hold on;
line([M_p,M_p], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_HFS_08192020(5,:),'FaceColor','none','LineWidth',1.2)
xlabel('M_p (\muA ms^{-1})')
ylabel('Number of Observations')
hold on;
line([M_p,M_p], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend(string(M_p)+' \muA ms^{-1}')

subplot(2,5,4)
hold on
line([p_p,p_p], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_HFS_08192020(6,:),'FaceColor','none','LineWidth',1.2)
xlabel('p_p (ms^{-1})')
ylabel('Number of Observations')
line([p_p,p_p], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend(string(p_p)+' ms^{-1}')


subplot(2,5,5)
hold on;
line([A_p,A_p], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_HFS_08192020(7,:),'FaceColor','none','LineWidth',1.2)
xlabel('A_p (\muA^2)')
ylabel('Number of Observations')
hold on;
line([A_p,A_p], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend(string(A_p)+' \muA^{2}')

subplot(2,5,8)
hold on;
line([M_d,M_d], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_HFS_08192020(8,:),'FaceColor','none','LineWidth',1.2)
xlabel('M_d (\muA ms^{-1})')
ylabel('Number of Observations')
hold on;
line([M_d,M_d], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend(string(M_d)+' \muA ms^{-1}')

subplot(2,5,9)
hold on;
line([p_d,p_d], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_HFS_08192020(9,:),'FaceColor','none','LineWidth',1.2)
xlabel('p_d (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([p_d,p_d], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend(string(p_d)+' ms^{-1}')

subplot(2,5,10)
hold on;
line([A_d,A_d], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_HFS_08192020(10,:),'FaceColor','none','LineWidth',1.2)
xlabel('A_d (\muA^2)')
ylabel('Number of Observations')
hold on;
line([A_d,A_d], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend(string(A_d)+' \muA^{2}')

%% Figure 23
clc

theta_LFS = mean(theta_LTD_08192020,2);

f = 281.3122;
gamma = theta_LFS(1);
eta = theta_LFS(2);
nu = theta_LFS(3);
g = theta_LFS(4);
M_p = theta_LFS(5);
p_p = theta_LFS(6);
A_p = theta_LFS(7);
M_d = theta_LFS(8);
p_d = theta_LFS(9);
A_d = theta_LFS(10);

std_gamma = std(theta_LTD_08192020(1,:));
std_eta = std(theta_LTD_08192020(2,:));
std_nu = std(theta_LTD_08192020(3,:));
std_g = std(theta_LTD_08192020(4,:));
std_M_p = std(theta_LTD_08192020(5,:));
std_p_p = std(theta_LTD_08192020(6,:));
std_A_p = std(theta_LTD_08192020(7,:));
std_M_d = std(theta_LTD_08192020(8,:));
std_p_d = std(theta_LTD_08192020(9,:));
std_A_d = std(theta_LTD_08192020(10,:));

figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
subplot(2,5,1)
hold on;
line([gamma,gamma], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_LTD_08192020(1,:),'FaceColor','none','LineWidth',1.2)
ylabel('Number of Observations')
xlabel('\gamma_C (ms^{-1})')
hold on;
line([gamma,gamma], [0, 160], 'LineWidth', 2, 'Color', 'r');
legend(string(gamma)+' ms^{-1}')
ylim([0,160])

subplot(2,5,2)
hold on;
line([eta,eta], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_LTD_08192020(2,:),'FaceColor','none','LineWidth',1.2)
xlabel('\eta (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([eta,eta], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend(string(eta)+' ms^{-1}')

subplot(2,5,6)
hold on;
line([nu,nu], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_LTD_08192020(3,:),'FaceColor','none','LineWidth',1.2)
xlabel('\nu (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([nu,nu], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend(string(nu)+' ms^{-1}')

subplot(2,5,7)
hold on;
line([g,g], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_LTD_08192020(4,:),'FaceColor','none','LineWidth',1.2)
xlabel('g (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([g,g], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend(string(g)+' ms^{-1}')

% figure(2)
subplot(2,5,3)
hold on;
line([M_p,M_p], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_LTD_08192020(5,:),'FaceColor','none','LineWidth',1.2)
xlabel('M_p (\muA ms^{-1})')
ylabel('Number of Observations')
hold on;
line([M_p,M_p], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend(string(M_p)+' \muA ms^{-1}')

subplot(2,5,4)
hold on
line([p_p,p_p], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_LTD_08192020(6,:),'FaceColor','none','LineWidth',1.2)
xlabel('p_p (ms^{-1})')
ylabel('Number of Observations')
line([p_p,p_p], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend(string(p_p)+' ms^{-1}')


subplot(2,5,5)
hold on;
line([A_p,A_p], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_LTD_08192020(7,:),'FaceColor','none','LineWidth',1.2)
xlabel('A_p (\muA^2)')
ylabel('Number of Observations')
hold on;
line([A_p,A_p], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend(string(A_p)+' \muA^{2}')

subplot(2,5,8)
hold on;
line([M_d,M_d], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_LTD_08192020(8,:),'FaceColor','none','LineWidth',1.2)
xlabel('M_d (\muA ms^{-1})')
ylabel('Number of Observations')
hold on;
line([M_d,M_d], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend(string(M_d)+' \muA ms^{-1}')

subplot(2,5,9)
hold on;
line([p_d,p_d], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_LTD_08192020(9,:),'FaceColor','none','LineWidth',1.2)
xlabel('p_d (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([p_d,p_d], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend(string(p_d)+' ms^{-1}')

subplot(2,5,10)
hold on;
line([A_d,A_d], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_LTD_08192020(10,:),'FaceColor','none','LineWidth',1.2)
xlabel('A_d (\muA^2)')
ylabel('Number of Observations')
hold on;
line([A_d,A_d], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend(string(A_d)+' \muA^{2}')


%% Figure 28

theta_DASOP = 10.^mean(theta_SOP_08182020,2);

clc

k1 = log10( theta_DASOP(1) );
k2 = log10( theta_DASOP(2) );
k3 = log10( theta_DASOP(3) );
k4 = log10( theta_DASOP(4) );
k5 = log10( theta_DASOP(5) );
k6 = log10( theta_DASOP(6) );
k7 = log10( theta_DASOP(7) );
k8 = log10( theta_DASOP(8) );

std_k1 = std(10.^(theta_SOP_08182020(1,:)));
std_k2 = std(10.^(theta_SOP_08182020(2,:)));
std_k3 = std(10.^(theta_SOP_08182020(3,:)));
std_k4 = std(10.^(theta_SOP_08182020(4,:)));
std_k5 = std(10.^(theta_SOP_08182020(5,:)));
std_k6 = std(10.^(theta_SOP_08182020(6,:)));
std_k7 = std(10.^(theta_SOP_08182020(7,:)));
std_k8 = std(10.^(theta_SOP_08182020(8,:)));



figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
subplot(2,4,1)
hold on;
line([k1,k1], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_SOP_08182020(1,:),'FaceColor','none','LineWidth',1.2)
xlabel('log(k_1) (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([k1,k1], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend('k_1 = '+string(10.^k1)+' ms^{-1}')


subplot(2,4,2)
hold on;
line([k2,k2], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_SOP_08182020(2,:),'FaceColor','none','LineWidth',1.2)
xlabel('log(k_2) (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([k2,k2], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend('k_2 = '+string(10.^k2)+' ms^{-1}')

subplot(2,4,3)
hold on;
line([k3,k3], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_SOP_08182020(3,:),'FaceColor','none','LineWidth',1.2)
xlabel('log(k_3) (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([k3,k3], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend('k_3 = '+string(10.^k3)+' ms^{-1}')

subplot(2,4,4)
hold on;
line([k4,k4], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_SOP_08182020(4,:),'FaceColor','none','LineWidth',1.2)
xlabel('log(k_4) (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([k4,k4], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend('k_4 = '+string(10.^k4)+' ms^{-1}')

subplot(2,4,5)
hold on;
line([k5,k5], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_SOP_08182020(5,:),'FaceColor','none','LineWidth',1.2)
xlabel('log(k_5) (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([k5,k5], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend('k_5 = '+string(10.^k5)+' ms^{-1}')

subplot(2,4,6)
hold on;
line([k6,k6], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_SOP_08182020(6,:),'FaceColor','none','LineWidth',1.2)
xlabel('log(k_6) (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([k6,k6], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend('k_6 = '+string(10.^k6)+' ms^{-1}')

subplot(2,4,7)
hold on;
line([k7,k7], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_SOP_08182020(7,:),'FaceColor','none','LineWidth',1.2)
xlabel('log(k_7) (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([k7,k7], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend('k_7 = '+string(10.^k7)+' ms^{-1}')

subplot(2,4,8)
hold on;
line([k8,k8], [0, 160], 'LineWidth', 2, 'Color', 'r');
histogram(theta_SOP_08182020(8,:),'FaceColor','none','LineWidth',1.2)
xlabel('log(k_8) (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([k8,k8], [0, 160], 'LineWidth', 2, 'Color', 'r');
ylim([0,160])
legend('k_8 = '+string(10.^k8)+' ms^{-1}')


%% Figure 29
theta_interaction = mean(theta_nov_24_2020,2);

clc

kHFS = theta_interaction(1); %k_sat
kmodstim = theta_interaction(2); %k_stim
knonlin = theta_interaction(3); %k_E
kLTD = theta_interaction(4) %k_late
knonlin_I  = theta_interaction(5); %k_I
tau_e = theta_interaction(6); %\tau_da
tau_d = theta_interaction(7); %\tua_stim
kDA = theta_interaction(8); %k_DA
kSOP = kHFS;

std_kHFS = std(theta_nov_24_2020(1,:));
std_kmodstim = std(theta_nov_24_2020(2,:));
std_knonlin = std(theta_nov_24_2020(3,:));
std_kLTD = std(theta_nov_24_2020(4,:));
std_knonlin_I  = std(theta_nov_24_2020(5,:));
std_tau_e = std(theta_nov_24_2020(6,:));
std_tau_d = std(theta_nov_24_2020(7,:));
std_kDA = std(theta_nov_24_2020(8,:));



figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
subplot(2,4,1)
hold on;
line([kHFS,kHFS], [0, 200], 'LineWidth', 2, 'Color', 'r');
histogram(theta_nov_24_2020(1,:),'FaceColor','none','LineWidth',1.2)
xlabel('k_{sat} (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([kHFS,kHFS], [0, 200], 'LineWidth', 2, 'Color', 'r');
ylim([0,200])
legend(string(kHFS)+' ms^{-1}')


subplot(2,4,2)
hold on;
line([kmodstim,kmodstim], [0, 200], 'LineWidth', 2, 'Color', 'r');
histogram(theta_nov_24_2020(2,:),'FaceColor','none','LineWidth',1.2)
xlabel('k_{stim} (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([kmodstim,kmodstim], [0, 200], 'LineWidth', 2, 'Color', 'r');
ylim([0,200])
legend(string(kmodstim)+' ms^{-1}')

subplot(2,4,3)
hold on;
line([knonlin,knonlin], [0, 200], 'LineWidth', 2, 'Color', 'r');
histogram(theta_nov_24_2020(3,:),'FaceColor','none','LineWidth',1.2)
xlabel('k_{E} (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([knonlin,knonlin], [0, 200], 'LineWidth', 2, 'Color', 'r');
ylim([0,200])
legend(string(knonlin)+' ms^{-1}')

subplot(2,4,4)
hold on;
line([kLTD,kLTD], [0, 200], 'LineWidth', 2, 'Color', 'r');
histogram(theta_nov_24_2020(4,:),'FaceColor','none','LineWidth',1.2)
xlabel('k_{late} (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([kLTD,kLTD], [0, 200], 'LineWidth', 2, 'Color', 'r');
ylim([0,200])
legend(string(kLTD)+'ms^{-1}')

subplot(2,4,5)
hold on;
line([knonlin_I,knonlin_I], [0, 200], 'LineWidth', 2, 'Color', 'r');
histogram(theta_nov_24_2020(5,:),'FaceColor','none','LineWidth',1.2)
xlabel('k_{I}')
ylabel('Number of Observations')
hold on;
line([knonlin_I,knonlin_I], [0, 200], 'LineWidth', 2, 'Color', 'r');
ylim([0,200])
legend(string(knonlin_I))

subplot(2,4,6)
hold on;
line([tau_e,tau_e], [0, 200], 'LineWidth', 2, 'Color', 'r');
histogram(theta_nov_24_2020(6,:),'FaceColor','none','LineWidth',1.2)
xlabel('\tau_{DR} (ms^{-1})')
ylabel('Number of Observations')
hold on;
line([tau_e,tau_e], [0, 200], 'LineWidth', 2, 'Color', 'r');
ylim([0,200])
legend(string(tau_e)+' ms^{-1}')

subplot(2,4,7)
hold on;
line([tau_d,tau_d], [0, 200], 'LineWidth', 2, 'Color', 'r');
histogram(theta_nov_24_2020(7,:),'FaceColor','none','LineWidth',1.2)
xlabel('\tau_{stim} (ms)')
ylabel('Number of Observations')
hold on;
line([tau_d,tau_d], [0, 200], 'LineWidth', 2, 'Color', 'r');
ylim([0,200])
legend(string(tau_d)+' ms')

subplot(2,4,8)
hold on;
line([kDA,kDA], [0, 200], 'LineWidth', 2, 'Color', 'r');
histogram(theta_nov_24_2020(8,:),'FaceColor','none','LineWidth',1.2)
xlabel('k_{DA} (ms)')
ylabel('Number of Observations')
hold on;
line([kDA,kDA], [0, 200], 'LineWidth', 2, 'Color', 'r');
ylim([0,200])
legend(string(kDA)+' ms')


%% Figure 20

load('theta_7best.mat')
load('theta_7best_NMDA.mat')

theta_f = mean(theta_7best,2);
tau1_ampa = theta_f(1); %ms
tau2_ampa = theta_f(2); %ms
gsbar = theta_f(3);

theta_f2 = mean(theta_7best_NMDA,2);
tau1_nmda  = theta_f2(1);
tau2_nmda  = theta_f2(2);
mew_nmda  = theta_f2(3);
gamma_nmda  = theta_f2(4);

figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
subplot(2,3,1)
hold on;
line([tau1_ampa,tau1_ampa], [0, 200], 'LineWidth', 2, 'Color', 'r');
histogram(theta_7best(1,:),'FaceColor','none','LineWidth',1.2)
% xlim([min(theta_7(1,:)) max(theta_7(1,:))])
xlabel('\tau_{1} (ms)')
ylabel('Number of Observations')
ylim([0,200])
hold on;
line([tau1_ampa,tau1_ampa], [0, 200], 'LineWidth', 2, 'Color', 'r');
legend(string(tau1_ampa)+' ms^{-1}')


subplot(2,3,2)
hold on;
line([tau2_ampa,tau2_ampa], [0, 200], 'LineWidth', 2, 'Color', 'r');
histogram(theta_7best(2,:),'FaceColor','none','LineWidth',1.2)
% xlim([min(theta_7(2,:)) max(theta_7(2,:))])
xlabel('\tau_{2} (ms)')
ylabel('Number of Observations')
ylim([0,200])
hold on;
line([tau2_ampa,tau2_ampa], [0, 200], 'LineWidth', 2, 'Color', 'r');
legend(string(tau2_ampa)+' ms^{-1}')

subplot(2,3,3)
hold on;
line([gsbar,gsbar], [0, 200], 'LineWidth', 2, 'Color', 'r');
histogram(theta_7best(3,:),'FaceColor','none','LineWidth',1.2)
% xlim([min(theta_7(3,:)) max(theta_7(3,:))])
xlabel('g_{s} (mS)')
ylabel('Number of Observations')
ylim([0,200])
hold on;
line([gsbar,gsbar], [0, 200], 'LineWidth', 2, 'Color', 'r');
legend(string(gsbar)+' mS')


% figure('DefaultAxesFontSize',14,'DefaultTextFontName','Calibri')
subplot(2,4,5)
hold on;
line([tau1_nmda,tau1_nmda], [0, 200], 'LineWidth', 2, 'Color', 'r');
histogram(theta_7best_NMDA(1,:),'FaceColor','none','LineWidth',1.2)
% xlim([min(theta_7(1,:)) max(theta_7(1,:))])
xlabel('\tau_{3} (ms)')
ylabel('Number of Observations')
ylim([0,200])
hold on;
line([tau1_nmda,tau1_nmda], [0, 200], 'LineWidth', 2, 'Color', 'r');
legend(string(tau1_nmda)+' ms^{-1}')

subplot(2,4,6)
hold on;
line([tau2_nmda,tau2_nmda], [0, 200], 'LineWidth', 2, 'Color', 'r');
histogram(theta_7best_NMDA(2,:),'FaceColor','none','LineWidth',1.2)
% xlim([min(theta_7(2,:)) max(theta_7(2,:))])
xlabel('\tau_{4} (ms)')
ylabel('Number of Observations')
ylim([0,200])
hold on;
line([tau2_nmda,tau2_nmda], [0, 200], 'LineWidth', 2, 'Color', 'r');
legend(string(tau2_nmda)+' ms^{-1}')

subplot(2,4,7)
hold on;
line([mew_nmda,mew_nmda], [0, 200], 'LineWidth', 2, 'Color', 'r');
histogram(theta_7best_NMDA(3,:),'FaceColor','none','LineWidth',1.2)
% xlim([min(theta_7(3,:)) max(theta_7(3,:))])
xlabel('\mu_s (\muM^{-1})')
ylabel('Number of Observations')
ylim([0,200])
hold on;
line([mew_nmda,mew_nmda], [0, 200], 'LineWidth', 2, 'Color', 'r');
legend(string(mew_nmda)+' \muM^{-1}')


subplot(2,4,8)
hold on;
line([gamma_nmda,gamma_nmda], [0, 200], 'LineWidth', 2, 'Color', 'r');
histogram(theta_7best_NMDA(4,:),'FaceColor','none','LineWidth',1.2)
% xlim([min(theta_7(4,:)) max(theta_7(4,:))])
xlabel('\gamma_s (mV^{-1})')
ylabel('Number of Observations')
ylim([0,200])
hold on;
line([gamma_nmda,gamma_nmda], [0, 200], 'LineWidth', 2, 'Color', 'r');
legend(string(gamma_nmda)+' mV^{-1}')

