
function PnMA= prob_nMA_old(xx,Px)

global LM LB PARAMS

% calculate how much we need to include in the EFOV
EFVOx= norminv(1 - PARAMS.I_FOV/2, 0, sqrt(Px(1,1)));
EFVOy= norminv(1 - PARAMS.I_FOV/2, 0, sqrt(Px(2,2)));
EFOV= sqrt(EFVOx^2 + EFVOy^2); 

% Get all visible landmarks, assuming no mis-extractions here
idf= get_visible_landmarks(xx,PARAMS.maxRange+EFOV, 0);

lm= LM(:,idf);
n_L= length(idf);

if n_L < 2, PnMA= 1; return, end;

% create table with yn
y2= NaN(n_L);
M= cell(n_L);
H= cell(n_L,1);
for t= 2:n_L
    [h_t,H{t}]= compute_lm_model(lm(:,t));
    
    for l= 1:(t-1)
        
        [h_l,H{l}]= compute_lm_model(lm(:,l));
        
        y= h_t - h_l;
        M{t,l}= (H{t} - H{l})*Px*(H{t} - H{l})';
        
        y2(t,l)= y'/M{t,l}*y;
    end    
end

y2star= zeros(n_L,1);
for t= 1:n_L
    
    [y2min,indMin]= min( [y2(t,1:t-1), NaN, y2(t+1:n_L,t)']  );
    sqrtM= sqrtm( M{max(t,indMin), min(t,indMin)} );
    
    Y= H{indMin}*Px*H{indMin}' + PARAMS.R; % assume same noise for all features
    
    lambda2= min( eig(sqrtM/Y*sqrtM) );
    
    % Lower bound weighted by M
    yStarM= sqrt(y2min) - sqrt( chi2inv(1 - PARAMS.I_y/n_L, PARAMS.m_F) );
    if yStarM > 0
        y2starM= (yStarM)^2;
    else 
        y2starM= 0;
    end
    
%     % Lower bound with the table
%     y2starM= interp1( LB.input, LB.table(n_L-1,:), y2min, 'linear','extrap');
    
    % Lower bound weighted by Y (non-centrality parameter)
    y2star(t)= y2starM*lambda2;
end


% Integrity risk
PnMA= 0;
for t= 1:n_L
    PnMA= PnMA + chi2cdf(0.25*y2star(t),5);
end
PnMA= (1 - PARAMS.I_y/n_L)*PnMA - (n_L - 1);

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
















































