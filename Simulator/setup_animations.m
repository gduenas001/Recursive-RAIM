function h= setup_animations()

global LM PARAMS


% Simulation Figure
% h.sim= figure('Units', 'pixels','position', [-1919, 1700, 1920, 220]); hold on; 
h.sim= figure('Units', 'pixels', 'position', [-1919+1919, 1700, 1920, 620]); hold on; 
xlabel('metres'), ylabel('metres')
% axis([0, 500, -20, 20]); axis equal
axis([-100, 40, -120, 20]); axis equal

% Initialize dynamic plots
h.xt= patch(0,0,'b'); % vehicle true
h.xv= patch(0,0,'r'); % vehicle estimate
h.xp= patch(0,0,'g'); % vehicle prediction estimate
h.pth= plot(0,0,'k.','markersize',2); % vehicle path estimate
h.obs= plot(0,0,'y'); % observations
h.vcov= plot(0,0,'r'); % vehicle covariance ellipses
h.fcov= plot(0,0,'r'); % feature covariance ellipses
h.path= plot(0,0,'r-'); % path of the vehicle
h.path_true= plot(0,0,'k-'); % True path

% Plot LMs
str= {};
for l= 1:size(LM,2)
    plot(LM(1,l),LM(2,l),'b.','markersize',12);
    str= [str, strcat('lm', num2str(l))];
    text(LM(1,l)-1,LM(2,l)-2,str(l),'FontSize',5,'color','b');
end


% The P(HMI)
h.figHMI= figure('Units', 'pixels','position', [-1919+1919 0 450 500]); hold on; grid on;
set(gca, 'fontsize', 10);
h.HMI= plot(0,0,'g-', 'linewidth',2);
h.HMI_CA= plot(0,0,'--b','linewidth',2);
legHMI= legend({'$P(HMI_k)$',...
    '$P(HMI_k | \neg F_K)$'},...
    'Interpreter','latex'); %, 'Location','northest');
set(legHMI, 'fontsize', 10);
set(gca,'Yscale','log');
xlim([1,PARAMS.numSteps]);
ylim([1e-12,1e-7]);
xlabel('Time epoch')


% Plot the P(nMA)_k and P(nMA)_K (at current time)
h.figPnMA= figure('Units', 'pixels','position', [-1400+1919 0 450 500]); hold on; grid on;
h.PMA_k= plot(0,0,'--b','linewidth',2);
h.PMA_K= plot(0,0,'-g','linewidth',2);
legend({'$P(MA_k | \neg MA_{K-1})$','$P(MA_K)$'}, 'Interpreter','latex');
xlim([1,PARAMS.numSteps]);


% Plot the epsilon error and cov envelope online
h.figEps= figure('Units', 'pixels','position', [-950+1919 0 450 500]); hold on; grid on;
set(gca, 'fontsize', 10);
h.eps= plot(0,0,'-b','linewidth',3);
h.cov= plot(0,0,'--g','linewidth',2);
legEps= legend({'$\hat{\epsilon}_k$','3-$\sigma$ Cov.'},'Interpreter','latex');
set(legEps, 'fontsize', 10);
xlim([1,PARAMS.numSteps]);
xlabel('Time epoch')
ylabel('meters');


% Plot the IAs online
h.figIA= figure('Units', 'pixels','position', [-500+1919 0 450 500]); hold on; grid on;
set(gca, 'fontsize', 10);
h.IA= plot(0,0,'-g','linewidth',2);
h.numAssoc= plot(0,0,'*k','linewidth',2);
legIA= legend({'$IA_k$', '# Associations'},'Interpreter','latex');
set(legIA, 'fontsize', 10);
xlim([1,PARAMS.numSteps]);
xlabel('Time epoch')
ylabel('');


% Plot the bias online
h.figBIAS= figure('Units', 'pixels','position', [-500+1919 0 450 500]); hold on; grid on;
set(gca, 'fontsize', 10);
h.BIAS= plot(0,0,'-g','linewidth',2);
legBIAS= legend({'$\epsilon_{\hat{H}}$'},'Interpreter','latex');
set(legBIAS, 'fontsize', 10);
xlim([1,PARAMS.numSteps]);
xlabel('Time epoch')
ylabel('');


% Plot the bias online
h.figRB= figure('Units', 'pixels','position', [-500+1919 0 450 500]); hold on; grid on;
set(gca, 'fontsize', 10);
h.q_RB= plot(0,0,'-g','linewidth',2);
h.T_RB= plot(0,0,'-r','linewidth',2);
legRB= legend({'$q_{RB}$', 'T_{RB}'},'Interpreter','latex');
set(legRB, 'fontsize', 10);
xlim([1,PARAMS.numSteps]);
xlabel('Time epoch')
ylabel('');



% % Plot the XX error and cov envelope online
% figErrorXX= figure('units','centimeters','pos', [60 0 15 10]); hold on; grid on;
% set(gca, 'fontsize', 10);
% h.errorXX= plot(0,0,'-b','linewidth',3);
% h.covXX= plot(0,0,'--g','linewidth',2);
% legErrorXX= legend({'$\hat{\epsilon}_k$ in X-axis','3-$\sigma$ Cov.'},'Interpreter','latex');
% set(legErrorXX, 'fontsize', 10);
% xlim([1,PARAMS.numSteps]);
% xlabel('Time epoch')
% ylabel('meters');




