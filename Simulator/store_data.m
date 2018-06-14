function store_data(step, P_HMI_noMA_UA, PnMA_k, xtrue, q_RB, T_RB)

global DATA XX BIAS PX

alpha= [-sin(XX(3)); cos(XX(3)); 0];
sig_hat= sqrt(alpha'*PX*alpha);

% Detector
DATA.q_RB(step)= q_RB;
DATA.T_RB(step)= T_RB;

% P(MA)
DATA.PnMA_k(step)= PnMA_k;
DATA.PnMA_K(step)= DATA.PnMA_K(step-1)*PnMA_k;

% HMI -- with UA
DATA.P_HMI_CA(step)= P_HMI_noMA_UA;
DATA.P_HMI(step)= 1 + ( P_HMI_noMA_UA - 1 )* DATA.PnMA_K(step);

% Error
DATA.epsXX(step,:)= abs(xtrue - XX)';
DATA.stdXX(step,:)= 3*sqrt(diag(PX)');

% Error state of interest
DATA.eps(step)= (alpha'* (xtrue - XX) );
DATA.stdEps(step)= 3*sig_hat;

% Worst case bias
DATA.bias(step,:)= BIAS;
DATA.bias_interest(step)= alpha'*BIAS;


% Path
DATA.path(step,:)= XX(1:2)';

