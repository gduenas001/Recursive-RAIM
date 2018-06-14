
function PnMA= prob_nMA(xx,Px)

global LM PARAMS
 
% % calculate how much we need to include in the EFOV
lambda_FOV= sqrt( max(eig(Px(1:2,1:2))) );
% EFOV= -lambda_FOV * norminv(PARAMS.I_FOV/4,0,1); % I think wrong!
EFOV= -sqrt(2) * lambda_FOV * norminv(PARAMS.I_FOV/2,0,1);
% EFOV= sqrt(2) * lambda_FOV * sqrt( chi2inv(1 - PARAMS.I_FOV,1) ); % same as previous




% Get all visible landmarks, assuming no mis-extractions here
idf= get_visible_landmarks(xx,PARAMS.maxRange+EFOV, 0);

lm= LM(:,idf);
n_L= length(idf);

if n_L < 2, PnMA= 1; return, end;



y2Star= ones(1,n_L).*inf;
for t= 1:n_L
    
    [h_t,H_t]= compute_lm_model(lm(:,t));
    Y= H_t * Px * H_t' + PARAMS.R;
    
    for l= 1:n_L
        if t == l, continue, end;
        
        [h_l,~]= compute_lm_model(lm(:,l));
        y= h_t - h_l;
        y2= y'/Y*y;
        
        if y2 < y2Star(t), y2Star(t)= y2; end;
    end
end


% Integrity risk
PnMA= 0; PnMA2= 0;
for t= 1:n_L
    PnMA= PnMA + chi2cdf(0.25*y2Star(t),5);
    PnMA2= PnMA2 + chi2cdf( ( sqrt(y2Star(t)) - sqrt(PARAMS.gate) )^2 , PARAMS.m_F);
end
PnMA= 1 - n_L - PARAMS.I_FOV + PnMA;

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
h= [d; atan2(dy,dx) - XX(3)];

% calculate H
H = [-dx/d, -dy/d,  0;
      dy/d2, -dx/d2, -1];



