
function [P_HMI_noMA_UA]= integrity_monitoring(remove_lm_ind)

global XX BIAS PX PARAMS

% Alpha using the prediction estimate
alpha= [-sin(XX(3)) ; cos(XX(3)) ; 0];

% Get all visible landmarks, assuming no mis-extractions here
lm_ind= get_visible_landmarks(XX, PARAMS.maxRange, 0);

% Remove the landmarks eliminated by LmS
lm_ind( ismember(lm_ind,remove_lm_ind) )= [];

% Number of landmarks to associate in the FV
n_L= length(lm_ind);
n= n_L * PARAMS.m_F;

% If no landmarks to associate --> return!
if n == 0
    l= PARAMS.alert_limit;
    sig_hat= alpha'*PX*alpha;
    P_HMI_noMA_UA= 1 - normcdf( (l + alpha'*BIAS)/sig_hat ) + normcdf( -(l + alpha'*BIAS)/sig_hat );
    return;
end

% RB Detector Threshold
T_RB= chi2inv(1 - PARAMS.C_REQ, n);

% Create model
H= zeros(n, 3);
Y= cell(n_L,1);
for l= 1:n_L
    ind= l*PARAMS.m_F-1; 
    [H(ind:ind+1,:), Y{l}]= lm_H_Y(lm_ind(l));
end
D= [H ; eye(3)];
R= kron(eye(n_L),PARAMS.R);
Delta= [R, zeros(n,3);
        zeros(3,n), PX];
invDelta= inv(Delta);
S= (D' * invDelta * D ) \ D' * invDelta; % S = I when no lms
DeltaI_DS= invDelta - invDelta * D * S;
 

% %% For the m states + some measurements worst-case fault
% 
% % Worst case direction
% E_f= [eye(n), zeros(n,3)];
% f_u= (E_f * (invDelta - invDelta * D * S)* E_f' ) \ E_f * S' * alpha;
% f_u= f_u / norm(f_u); % normalize -- I don't know if necessary
% 
% % Worst case fault norm
% PX_hat= S*Delta*S';
% sig_hat= sqrt(alpha'*PX_hat*alpha);
% 
% % % Worst case fault norm for each UA
% % f_norm_max= -1;
% % for l= 1:n_L
% %     ind= l*PARAMS.m_F - 1;
% %     f_norm_max_l= individual_worst_fault_norm(H(ind:ind+1,:),Y{l},f_u(ind:ind+1));
% %     if f_norm_max_l > f_norm_max
% %         f_norm_max= f_norm_max_l;
% %     end
% % end
% 
% % Minimize function
% P_HMI_noMA_UA_fn= @(f_norm) (-1)*...
%     ( 1 - cdf('normal', PARAMS.alert_limit, -alpha' * S * [f_u*f_norm; 0; -alpha'*BIAS; 0], sig_hat) + ...
%     cdf('normal', -PARAMS.alert_limit, -alpha' * S * [f_u*f_norm; 0; -alpha'*BIAS; 0], sig_hat) ) * ...
%     cdf('Noncentral Chi-square', T_RB, n, [f_u*f_norm; 0; -alpha'*BIAS; 0]' * DeltaI_DS * [f_u*f_norm; 0; -BIAS(2); 0]);
% 
% f_norm= fminbnd(P_HMI_noMA_UA_fn, -4, 4);
% 
% % If P(HMI) = 0 for all fault norms -> keep f_norm at 0
% % if abs(P_HMI_noMA_UA_fn(f_norm) - P_HMI_noMA_UA_fn(0)) < 1e-20
% %     f_norm= 0;
% % end
% 
% 
% % HMI with no MA
% P_HMI_noMA_UA= -P_HMI_noMA_UA_fn(f_norm);
% 
% % Update bias
% BIAS= - S * [f_u*f_norm; 0; -BIAS(2); 0];


%% Assume only previous faults affect
E_feps= [zeros(3,n), eye(3,3)];
feps_u= (E_feps * (DeltaI_DS)* E_feps' ) \ E_feps * S' * alpha;
feps_u= feps_u / norm(feps_u); % normalize -- I don't know if necessary
f_l_u= [zeros(n,1); feps_u];

% Worst case fault norms
PX_hat= S*Delta*S';
sig_hat= sqrt(alpha'*PX_hat*alpha);
fhat_eps= -alpha' * S * f_l_u;
lambda2= f_l_u' * DeltaI_DS * f_l_u;


P_HMI_noMA_UA_fn= @(f_norm) (-1)* ( 1 - cdf('normal', PARAMS.alert_limit, fhat_eps*f_norm, sig_hat) + ...
    cdf('normal', -PARAMS.alert_limit, fhat_eps*f_norm, sig_hat) ) * ...
    cdf('Noncentral Chi-square', T_RB, n, lambda2*f_norm^2);

[f_norm,P_HMI_noMA_UA]= fminbnd(P_HMI_noMA_UA_fn, 0, 4);

P_HMI_noMA_UA= -P_HMI_noMA_UA;
BIAS= - S * f_l_u*f_norm;



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% function [P_HMI_noMA_UA]= P_HMI_noMA_UA_fn(f_norm, f_u, n_L, sig_hat, S, alpha, DeltaI_DS, T_RB)
% 
% global PARAMS BIAS
% 
% n= n_L*PARAMS.m_F;
% f= zeros(n,1);
% for i= 1:n_L
%     ind= (i-1)*PARAMS.m_F + 1;
%     inds= ind:(ind+1);
%     
%     f(inds)= f_norm(i)*f_u(inds);
% end
% f_l= [f; -BIAS];
% 
% % Fault mean
% f_eps= -alpha' * S * f_l;
% 
% % Detector non-centrality paramter
% lambda2= f_l' * DeltaI_DS * f_l;
% 
% P_HMI_noMA_UA= ( 1 - cdf('normal', PARAMS.alert_limit, f_eps, sig_hat) + ...
%     cdf('normal', -PARAMS.alert_limit, f_eps, sig_hat) ) * ...
%     cdf('Noncentral Chi-square', T_RB, n, lambda2);








% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [H,Y]= lm_H_Y(lm_id)
% return normalised innovation squared (ie, Mahalanobis distance)

global XX PX LM PARAMS

% auxiliary values
dx= LM(1,lm_id) - XX(1); 
dy= LM(2,lm_id) - XX(2);
d2= dx^2 + dy^2;
d= sqrt(d2);

xd= dx/d;
yd= dy/d;
xd2= dx/d2;
yd2= dy/d2;

% calculate H
H = [-xd -yd 0; 
      yd2 -xd2 -1];

Y= H*PX*H' + PARAMS.R;  

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [f_max]= individual_worst_fault_norm(H,Y,f_u)

global BIAS PARAMS

num= sqrt(PARAMS.T2) + (H*BIAS)'/Y*H*BIAS;
den= f_u'/Y*f_u;
f_max= num / den;








