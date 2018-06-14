
function post_processing_and_plots(step)

global DATA SWITCH

% clean extra allocated memory
DATA.PCA(step+1:end)= [];
DATA.PCAt(step+1:end)= [];
DATA.realPCA(step+1:end)= [];
DATA.calcPCA(step+1:end)= [];
DATA.errorXX(step+1:end,:)= [];
DATA.stdXX(step+1:end,:)= [];
DATA.PHMI(step+1:end)= [];
DATA.numLmAssoc(step+1:end)= [];
DATA.numMsmts(step+1:end)= [];


% plots - P(CA) VS time
figure; hold on; grid on;
title('Instant P(CA)');
xlabel('Time'); ylabel('P(CA)');
plot(DATA.PCA,'-b','linewidth',3);
% plot(DATA.PCAt,'or');
% idx= find(DATA.PCAt == 0);
% for i= 1:length(idx)
%     line([idx(i),idx(i)],[0,1],'color','red');
% end
% axis([0,step,0,1]);


% plots - Error Vs Covariance - X & Y
figure; hold on; grid on
plot(DATA.errorXX(:,1),'b-');
plot(DATA.errorXX(:,2),'r-');
plot(DATA.stdXX(:,1),'b--','linewidth',5);
plot(DATA.stdXX(:,2),'r--','linewidth',5);
% for i= 1:length(idx)
%     line([idx(i),idx(i)],[0,max(DATA.stdXX(:))],'color','black');
% end
xlabel('time')
ylabel('m')
legend('error X','error Y','covariance X','covariance Y','location','northeast');

% % plots - Error Vs Covariance - phi
% figure; hold on; grid on
% plot(errorXX(:,3),'g-');
% plot(stdXX(:,3),'g--','linewidth',5);
% for i= 1:length(idx)
%     line([idx(i),idx(i)],[0,max(stdXX(:,3))],'color','black');
% end
% xlabel('time')
% ylabel('rad')
% legend('error phi','covariance phi','location','southeast');

% % plots - realPCA Vs calcPCA
% subplot(1,2,2); hold on; grid on
% title('Averaged P(CA)');
% xlabel('Time'); ylabel('Averaged P(CA)');
% % plot(DATA.realPCA,'-r','linewidth',2);
% plot(DATA.calcPCA,'b-','linewidth',3);
% plot(DATA.calcPCA_MJ,'g--','linewidth',3);
% % legend('Real P(CA)','Average calculated P(CA)','Average calculated P(CA) - MJ','location','southeast');
% axis([0,step,0,1]);
% 


% find(errorXX(:,1) > stdXX(:,1))
% find(errorXX(:,2) > stdXX(:,2))
% find(errorXX(:,3) > stdXX(:,3))


if SWITCH.association == 0 % known DA
    
    gamma_plot1= zeros(2,step);
    gamma_plot2= zeros(2,step);
    stdngamma_plot1= zeros(2,step);
    stdngamma_plot2= zeros(2,step);
    
    for i= 1:step
        gamma_plot1(1:2,i)= DATA.gamma{i}(1:2);
        gamma_plot2(1:2,i)= DATA.gamma{i}(3:4);
        stdngamma_plot1(1:2,i)= DATA.stdngamma{i}(1:2);
        stdngamma_plot2(1:2,i)= DATA.stdngamma{i}(3:4);
    end
    
    figure; title('Innovation vectors for each landmark');
    subplot(1,2,1); hold on; grid on;
    plot(gamma_plot1(1,:),gamma_plot1(2,:),'o');
    xlabel('range innovation');
    ylabel('angle innovation')
    axis equal
    
    subplot(1,2,2); hold on; grid on;
    plot(gamma_plot2(1,:),gamma_plot2(2,:),'o');
    xlabel('range innovation');
    ylabel('angle innovation')
    axis equal
    
    figure; title('Standard Innovation vectors for each landmark');
    subplot(1,2,1); hold on; grid on;
    plot(stdngamma_plot1(1,:),stdngamma_plot1(2,:),'o');
    xlabel('range innovation');
    ylabel('angle innovation')
    axis equal
    
    subplot(1,2,2); hold on; grid on;
    plot(stdngamma_plot2(1,:),stdngamma_plot2(2,:),'o');
    xlabel('range innovation');
    ylabel('angle innovation')
    axis equal
end

if SWITCH.association == 3 || SWITCH.association == 4
    figure; hold on; grid on;
    title('P(HMI)')
    plot(DATA.PHMI,'b');
    xlabel('Time');
%     axis([0,step,0,1]);
    
    figure; hold on; grid on;
    title('Number of msmts/lms associated Vs number of msmts');
    plot(DATA.numLmAssoc,'*r');
    plot(DATA.numMsmts, 'g-');
    xlabel('Time');
end

























