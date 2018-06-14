
function [PnMA,remove_lm_index,EFV,n_L_LS]= prob_nMA_LS(xx,Px)

global LM PARAMS
 
% % calculate how much we need to include in the EFOV
lambda_FV= sqrt( max(eig(Px(1:2,1:2))) );
EFV= -sqrt(2) * lambda_FV * norminv(PARAMS.I_FOV/2,0,1);
% EFOV= sqrt(2) * lambda_FOV * sqrt( chi2inv(1 - PARAMS.I_FOV,1) ); % same as previous


% Get all visible landmarks, assuming no mis-extractions here
idf= get_visible_landmarks(xx,PARAMS.maxRange+EFV, 0);

lm= LM(:,idf);
n_L= length(idf);

% Less than 2 lms -> no possible MA
if n_L < 2
    PnMA= 1;
    remove_lm_index= [];
    n_L_LS= n_L;
    return
end 

% LS threshold
T2_y= 4* chi2inv( 1 + (PARAMS.I_FOV-PARAMS.P_IA_max)/(n_L), PARAMS.m_F);
% T2_y= ( sqrt( chi2inv( 1 + (PARAMS.I_FOV-PARAMS.P_IA_max)/(n_L), PARAMS.m_F) ) + sqrt(PARAMS.T2) )^2;

remove_lm_index= [];
y2Star= ones(1,n_L).*inf;
for t= 1:n_L
    [h_t,H_t]= compute_lm_model(lm(:,t));
    Y= H_t * Px * H_t' + PARAMS.R;
    
    for l= t+1:n_L
        
        [h_l,~]= compute_lm_model(lm(:,l));
        y= h_t - h_l; y(2)= pi_to_pi(y(2));
        y2= y'/Y*y;
        
        if y2 < T2_y
            remove_lm_index= [remove_lm_index, idf(t), idf(l)];
        elseif y2 < y2Star(t)
            y2Star(t)= y2;
        end
    end
end

remove_lm_index= unique(remove_lm_index);
n_L_LS= n_L - length(remove_lm_index);

% Integrity risk
PnMA= 0;
for t= 1:n_L
    if ~ismember(idf(t),remove_lm_index)
        PnMA= PnMA + chi2cdf(0.25*y2Star(t), PARAMS.m_F);
%         PnMA= PnMA + chi2cdf( ( sqrt(y2Star(t)) - sqrt(PARAMS.T2) )^2 , PARAMS.m_F);
%         PnMA= PnMA + chi2cdf(...
%             ( 0.5*sqrt(y2Star(t)) - (BIAS'*H_t{t}') / Y{t} * (H_t{t}*BIAS) )^2 ...
%             ,PARAMS.m_F);
    end
end
PnMA= 1 - n_L_LS - PARAMS.I_FOV + PnMA;

% Eliminate negative values;
PnMA= max(PnMA,0);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h,H]= compute_lm_model(lm)

global XX

dx= lm(1) - XX(1);
dy= lm(2) - XX(2);
d2= dx^2 + dy^2;
d= sqrt(d2);

% calculate h
h= [d; 
    pi_to_pi( atan2(dy,dx) - XX(3) )];

% calculate H
H = [-dx/d, -dy/d,  0;
      dy/d2, -dx/d2, -1];



