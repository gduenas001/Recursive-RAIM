
clear; clc; close all;


load('DATA.mat');

numSteps= 1300;


DATA.P_HMI_CA( DATA.P_HMI_CA == 0 )= 1e-20;


%% Integrity risk
figure; hold on; grid on;
plot(1:numSteps ,DATA.P_HMI(1:numSteps), 'g', 'linewidth',2);
plot(1:numSteps , max(eps,DATA.P_HMI_CA(1:numSteps)), 'b', 'linewidth',2);
xlim([0,numSteps]);
ylim([eps,1]);
xlabel('Time epoch')
ylabel('Probability')
set(gca, 'fontsize', 10);
set(gca,'Yscale','log');
ax= gca;
ax.YTick= [1e-14, 10e-13, 10e-11, 10e-9, 10e-7, 10e-5, 10e-3, 10e-1];
ax.XTick= [0:200:1200];

legdHMI= legend({'$\breve{P}(HMI_k)$',...
    '$P(HMI_k ~|~ \neg IA_K)$'},...
    'Interpreter','latex', 'Location','northwest');
legdHMI.FontSize= 12;

%% Error
figure; hold on; grid on;
plot(1:numSteps ,DATA.eps(1:numSteps), 'r', 'linewidth',2);
plot(1:numSteps ,DATA.stdEps(1:numSteps), '-b', 'linewidth',2);
plot(1:numSteps ,-DATA.stdEps(1:numSteps), '-b', 'linewidth',2);
xlim([0,numSteps]);
ylim([-5.7,1.5]);
xlabel('Time epoch')
ylabel('meters')
set(gca, 'fontsize', 10);

legdEps= legend({'$\hat{\epsilon}_k$ Estimate Error',...
    '$3 \hat{\sigma}_k$ Covariance Envelope'},...
    'Interpreter','latex', 'Location','southwest');
legdEps.FontSize= 12;
































