function [P_HMI_worst,F_mag,Fault_slope,H_M,L_M,L_p_M,L_pp_M,Y]= integrity_monitoring_fault...
    (Px, Phi_M ,H_M ,L_M ,L_p_M ,L_pp_M ,A_k ,Y ,T , P_CA)

% Px : states prediction covarience matrix
% Phi_M : state tranision matrix over the horizon, including the current (concatenated)
% H_M : observation matrix over the horizon, excluding the current (concatenated)
% L_M : Kalman gain over the horizon, excluding the current (concatenated)
% L_p_M : L_prime over the horizon, excluding the current (concatenated)
% L_pp_M : L_prime_prime over the horizon, excluding the current (concatenated)
% A_k : () current
% Y : Innovation vector covarience matrix during the horizon
% P_CA : The previous probability of correct asscoiations
% T : Gate threshold

global LM PARAMS
M= PARAMS.Preceding_Horizon;
m_F= PARAMS.m_F;

% % calculate how much we need to include in the EFOV
%lambda_FV= sqrt( max(eig(Px(1:2,1:2))) );
%EFV= -sqrt(2) * lambda_FV * norminv(PARAMS.I_FOV/2,0,1);
% EFOV= sqrt(2) * lambda_FOV * sqrt( chi2inv(1 - PARAMS.I_FOV,1) ); % same as previous


% Get all visible landmarks, assuming no mis-extractions here
%idf= get_visible_landmarks(xx,PARAMS.maxRange+EFV, 0);

%lm= LM(:,idf);
%n_L= length(idf);

% Assuming that LIDAR range is infinte, and there are no misextractions
n_L= size(LM,2);

H_k=zeros(n_L*m_F, PARAMS.m); % Observation matrix at the current time step
h_k=zeros(n_L*m_F, 1); % Expected measurement
Y_k=zeros(n_L*m_F, m_F); % Current IIV covarience (concatenated)

for t= 1:n_L
    [h_t,H_t]= compute_lm_model(LM(:,t));
    H_k(((t-1)*m_F)+1:t*m_F,:)= H_t;
    h_k(((t-1)*m_F)+1:t*m_F,:)= h_t;
    Y_k(((t-1)*m_F)+1:t*m_F,:)= H_t * Px * H_t' + PARAMS.R;
end

min_y_norm= zeros(n_L,1); % the minimum weighted norm of y for every landmark (concatenated)


for t= 1:n_L
    dummy_variable= inf;
    h_t= h_k(((t-1)*m_F)+1:t*m_F,:);
    Y_t= Y_k(((t-1)*m_F)+1:t*m_F,:);
    % Find the smallest separation for each landmark
    for j=1:n_L
        if j == t, continue; 
        else
            h_j= h_k(((j-1)*m_F)+1:j*m_F,:);
            y_j= h_t - h_j;
            Norm_y_j= sqrt(y_j'*Y_t*y_j);
            dummy_variable= min(Norm_y_j,dummy_variable);
        end
    end
    min_y_norm(t)= dummy_variable;
end


% % Less than 2 lms -> no possible MA
% if n_L < 2
%     PnMA= 1;
%     remove_lm_index= [];
%     n_L_LS= n_L;
%     return
% end 

% % LS threshold
% T2_y= 4* chi2inv( 1 + (PARAMS.I_FOV-PARAMS.P_IA_max)/(n_L), m_F);
% % T2_y= ( sqrt( chi2inv( 1 + (PARAMS.I_FOV-PARAMS.P_IA_max)/(n_L), m_F) ) + sqrt(PARAMS.T2) )^2;
% 
% remove_lm_index= [];
% y2Star= ones(1,n_L).*inf;
% H_k=[];
% for t= 1:n_L
%     [h_t,H_t]= compute_lm_model(lm(:,t));
%     Y= H_t * Px * H_t' + PARAMS.R;
%     
%     for l= t+1:n_L
%         
%         [h_l,~]= compute_lm_model(lm(:,l));
%         y= h_t - h_l; y(2)= pi_to_pi(y(2));
%         y2= y'/Y*y;
%         
%         if y2 < T2_y
%             remove_lm_index= [remove_lm_index, idf(t), idf(l)];
%         elseif y2 < y2Star(t)
%             y2Star(t)= y2;
%             H_k=[H_k;H_t];
%         end
%     end
% end


V= kron(eye(n_L),PARAMS.R); % current measurements covarience matrix
Y_k= H_k * Px * H_k' + V; % Current IV covarience matrix

% Update the innovation vector covarience matrix for the new PH
Y(m_F*n_L+1:end,m_F*n_L+1:end)= Y(1:end-m_F*n_L,1:end-m_F*n_L);
Y(1:m_F*n_L,1:m_F*n_L)= Y_k;

Lk= Px * H_k / Y_k;
P_Hat= Px - Lk*H_k*Px;

Lk_p= eye(PARAMS.m) - Lk*H_k; % current ()
Lk_pp= Lk_p * Phi_M(1:PARAMS.m,1:PARAMS.m); % current ()

% Update the horizon of L matrix (Kalman Gain)
L_M= [Lk;L_M];
L_p_M= [Lk_p;L_p_M];
L_pp_M= [Lk_pp;L_pp_M];

% Update the horizon of H matrix
H_M= [H_k;H_M];

%remove_lm_index= unique(remove_lm_index);
%n_L_LS= n_L - length(remove_lm_index);

if isempty(A_k)
    
    % Create A_k^M for the first time (later it will be done recursively)
    A_k= inf* ones( PARAMS.m, n_L*m_F*(M+1) + PARAMS.m );
    A_k(: , 1:n_L*m_F)= Lk;
    for i= 1:M
        if i == 1
            Dummy_Variable= L_pp_M(1:i*PARAMS.m,:);
        else
            Dummy_Variable= Dummy_Variable*L_pp_M( (i-1)*PARAMS.m+1:i*PARAMS.m, : );
        end
        A_k(:, n_L*m_F*i + 1:2*n_L*m_F*i)=...
            Dummy_Variable * L_M( i*n_L*m_F+1 : 2*i*n_L*m_F, : );
    end
    A_k( :,n_L*m_F*(M+1)+1 : 2*n_L*m_F*(M+1) )=...
        Dummy_Variable * L_pp_M( M*PARAMS.m + 1 : (M+1)*PARAMS.m, : );
    
    % current()
    %B_k=[eye(n_L*m_F),-H_k*G(1:PARAMS.m,1:PARAMS.m)*Lk_pp\A_k(:,(n_L*m_F)+1:end)];
    
    % augmented B (map faults during horizon to innovations during horizon)
    B_bar= inf* ones( n_L*m_F*(M+1) : n_L*m_F*(M+1)+PARAMS.m );
    Dummy_variable= A_k;
    Dummy_variable= Lk_pp \ Dummy_variable( : , n_L*m_F + 1:end ); % A_(k-1)^(M-1)
    B_bar(1:n_L*m_F,:)= [eye(n_L*m_F), -H_k*Phi_M(1:PARAMS.m,:)*Dummy_variable]; % B_k^M
    % Recursive computation of B
    for i= 1:M
        Dummy_variable= L_pp_M( i*PARAMS.m+1:i*PARAMS.m + PARAMS.m , :)...
            \ Dummy_variable(:,(n_L*m_F)+1:end);
        B_bar(i*n_L*m_F+1:2*i*n_L*m_F,:)= ...
            [zeros(n_L*m_F,n_L*m_F*i),eye(n_L*m_F),-H_M(n_L*m_F*i+1:2*n_L*m_F*i,:)*Phi_M(PARAMS.m*i+1:2*PARAMS.m*i,:)*Dummy_variable];
    end
    M_k= B_bar' / Y * B_bar;

else % The previous A_(k-1) is an input to the function -> A_k is computed recursively
    
    % Calculating A_k recursively
    A_k= [Lk, Lk_pp*A_k(:,1:n_L*m_F*M), Lk_pp*A_k(:,end-PARAMS.m+1:end) / L_pp_M(:,end-m+1:end)];
    
    % current()
    %B_k=[eye(n_L*m_F),-H_k*G(1:PARAMS.m,1:PARAMS.m)*Lk_pp\A_k(:,(n_L*m_F)+1:end)];
    
    % B_bar - SAME AS BEFORE (CREATE A FUNCTION)
    B_bar= inf* ones( n_L*m_F*(M+1) : n_L*m_F*(M+1)+PARAMS.m );
    Dummy_variable= A_k;
    Dummy_variable= Lk_pp \ Dummy_variable( : , n_L*m_F + 1:end ); % A_(k-1)^(M-1)
    B_bar(1:n_L*m_F,:)= [eye(n_L*m_F), -H_k*Phi_M(1:PARAMS.m,:)*Dummy_variable]; % B_k^M
    % Recursive computation of B
    for i= 1:M
        Dummy_variable= L_pp_M( i*PARAMS.m+1:i*PARAMS.m + PARAMS.m , :)...
            \ Dummy_variable(:,(n_L*m_F)+1:end);
        B_bar(i*n_L*m_F+1:2*i*n_L*m_F,:)= ...
            [zeros(n_L*m_F,n_L*m_F*i),eye(n_L*m_F),-H_M(n_L*m_F*i+1:2*n_L*m_F*i,:)*Phi_M(PARAMS.m*i+1:2*PARAMS.m*i,:)*Dummy_variable];
    end
    M_k= B_bar' / Y * B_bar;
    
    % Removing the last element (matrix) to fix the horizon
    H_M= H_M(1:end-(n_L*m_F),:);
    L_M= L_M(1:end-m,:);
    L_p_M= L_p_M(1:end-m,:);
    L_pp_M= L_pp_M(1:end-m,:);
end




% Initialization of fault slope for every hypothesis (concatenated)
% Hypotheses include both the cases where there are no previous faults in 
% the states, and the cases where there are previous faults (the fault free
% case for both initial states, and horizon measurements are not included)
% Organization:
% (only previous state faults ;;;; 
% both previous state faults, and the whole combination of single measurement faults ;;;; 
% previous states are fault-free and the whole combination of single measurement faults)

Fault_slope= inf* ones( (n_L*m_F*(M+1)+PARAMS.m)*(n_L*(M+1))*2+1, 1 );

% E matrix for only previous state faults
E= zeros( m, (n_L*m_F*(M+1)+PARAMS.m) );
E(:, end-m+1:end)= eye(PARAMS.m);

% Worst-case fault slope for only previous state faults
Fault_slope( 1 : n_L*m_F*(M+1)+PARAMS.m )= E'/(E*M_k*E')*E*A_k'*alpha;

% Loop over hypotheses with faulted measurements in the PH (only 1 fault)
for i= 1:2*n_L*(M+1)
    
    if i <= n_L*(M+1) % Hypotheses with previous faults
        
        % E matrix for both previous state faults,and the whole combination
        % of single measurement faults
        E= zeros( PARAMS.m + m_F , n_L*m_F*(M+1) + PARAMS.m );
        E( end-PARAMS.m+1 : end,end-PARAMS.m+1:end )= eye(PARAMS.m);
        E( 1:m_F , (i-1)*m_F+1 : i*m_F )= eye(m_F);
        
        % Worst-case fault slope for both previous state faults,and the 
        % whole combination of single measurement faults
        Fault_slope( (n_L*m_F*(M+1)+PARAMS.m)*i+1 : 2*i*(n_L*m_F*(M+1)+PARAMS.m) )=...
            E'/(E*M_k*E')*E*A_k'*alpha;
    
    else % Hypotheses without previous faults
        
        % E matrix for previous states are fault-free,and the whole 
        % combination of single measurement faults
        E= zeros(m_F , n_L*m_F*(M+1) + PARAMS.m);
        E( : , (i-1-(M+1))*m_F+1:(i-(M+1))*m_F )= eye(m_F);
        
        % Worst-case fault slope for previous states are fault-free,and the
        % whole combination of single measurement faults
        Fault_slope( (n_L*m_F*(M+1)+PARAMS.m)*i+1:2*i*(n_L*m_F*(M+1)+PARAMS.m) )=...
            E'/(E*M_k*E')*E*A_k'*alpha;
    
    end
end

% Assuming that fault free case for the calculation of probability of
% correct associations
if isempty(P_CA) == 1, P_CA= 1; end
Dummy_Variable_2= 1-n_L;
for i= 1:n_L
    Dummy_Variable_2 = Dummy_Variable_2 +...
       cdf( 'Central Chi-square' , max( (0.5*min_y_hat(i))^2 , (min_y_hat(i)-T)^2 ) );
end
P_CA_k= Dummy_Variable_2 * P_CA;

% if (isempty(P_CA) == 1)
%     P_CA=ones(2*(n_L*(M+1))+2,1);
% end
% P_CA_k=zeros(2*(n_L*(M+1))+2,1);
% 
% Dummy_Variable_2=1-n_L;
% for i=1:n_L
%     Dummy_Variable_2 = Dummy_Variable_2 + cdf( 'Central Chi-square' , max( (0.5*min_y_hat(i))^2 , (min_y_hat(i)-T)^2 ) );
% end
% 
% P_CA_k(1) = Dummy_Variable_2 * P_CA(1);
% 
% for j=1:2*(n_L*(M+1))+1
%     
%     Dummy_Variable_2=1-n_L;
%     Fault=Fault_mag * Fault_slope((n_L*m_F*(M+1)+PARAMS.m)*(j-1)+1:(n_L*m_F*(M+1)+PARAMS.m)*j);
%     
%     for i=1:n_L
%         Y_i=Y_k(((i-1)*m_F)+1:i*m_F,:);
%         Norm_H_ti_F_epsi= sqrt( Fault' * [ zeros(m_F*n_L) , B_k(:,m_F*n_L+1 : end)]' * [ zeros(m_F,(i-1)*m_F) , eye(m_F) , zeros(m_F,n_L*m_F-i*m_F) ]'* Y_i * [ zeros(m_F,(i-1)*m_F) , eye(m_F) , zeros(m_F,n_L*m_F-i*m_F) ] * [ zeros(m_F*n_L) , B_k(:,m_F*n_L+1 : end)] * Fault );
%         Dummy_Variable_2 = Dummy_Variable_2 + cdf( 'Central Chi-square' , max( (0.5*min_y_hat(i) - Norm_H_ti_F_epsi)^2 , (min_y_hat(i) - Norm_H_ti_F_epsi - T)^2 ));
%     end
%     P_CA_k(j+1) = Dummy_Variable_2 * P_CA(j+1);
% end

syms Fault_mag

% Detector threshold including the PH
T_RB= chi2inv(1 - PARAMS.C_REQ, n_L*m_F*(M+1));

% First term in integrity risk equation for every hypothesis (concatenated)
P_HMI_H= inf* ones(2*(n_L*(M+1))+2,1); 
P_HMI_H(1) = (1-cdf('Normal',PARAMS.alert_limit,0,alpha*P_Hat*alpha')+cdf('Normal',-PARAMS.alert_limit,0,alpha*P_Hat*alpha')) * cdf('Central Chi-square',T_RB,m_F*n_L*(M+1));
% Loop over every hypothesis
for i= 1:2*(n_L*(M+1))+1    
    mu=alpha*A_k*Fault_mag*Fault_slope((n_L*m_F*(M+1)+PARAMS.m)*(i-1)+1:(n_L*m_F*(M+1)+PARAMS.m)*i);
    non_centrality=(Fault_mag^2)*(Fault_slope((n_L*m_F*(M+1)+PARAMS.m)*(i-1)+1:(n_L*m_F*(M+1)+PARAMS.m)*i))'*B_bar'/Y*B_bar*(Fault_slope((n_L*m_F*(M+1)+PARAMS.m)*(i-1)+1:(n_L*m_F*(M+1)+PARAMS.m)*i));
    P_HMI_H(i+1) = (1-cdf('Normal',PARAMS.alert_limit,mu,alpha*P_Hat*alpha')+cdf('Normal',-PARAMS.alert_limit,mu,alpha*P_Hat*alpha')) * cdf('Noncentral Chi-square',T_RB,m_F*n_L*(M+1),non_centrality);
end

% Initialization
P_HMI= 0;

% Compute the P(HMI) given all the terms
for i= 1:2*(n_L*(M+1))+2
    P_HMI=P_HMI+(1+(P_HMI_H(i)-1)*P_CA_k)*PARAMS.P_UA;
end


[F_mag, P_HMI_worst]= fminbnd(@(Fault_mag) -P_HMI, -10, 10);

P_HMI_worst=-P_HMI_worst;

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
