

function  generate_models(z,dz,Nz,R_2x2,V_FOV, LAMBDA, P_D)

global T h hlm H Hlm gamma ngamma Y psi beta phi Const
global XX PX PARAMS

% compute the constant C
Const= (2*pi)^(dz/2) * P_D^(-1) * LAMBDA * V_FOV;
log_C2= log(Const^2);

beta= zeros(psi,1);
phi= zeros(1,psi);
H= cell(psi,1);
h= cell(psi,1);
gamma= cell(psi,1);
ngamma= zeros(psi,1);
Y= cell(psi,1);
for j= 1:psi
    
    % current possible association
    T_j= T(j,:);
    phi_j= T_j(end);
    
    % If all outliers
    if phi_j == Nz, beta(j)= -phi_j*log_C2; continue, end;
    
    % If there is an outlier eliminate the measurement associated with it
    zj= z(:, T_j(1:end-1) ~= 0 );
    
    % make R a block diagonal with higher dimension
    R= kron(eye(Nz-phi_j),R_2x2);
    
    % create the models for each association
    for i= 1:Nz
        if T_j(i) ~= 0
            H{j}= [H{j};Hlm{T_j(i)}];
            h{j}= [h{j};hlm{T_j(i)}];
        end
    end
    gamma{j}= zj(:) - h{j};
    Y{j}= H{j}*PX*H{j}' + R;
    ngamma(j)= gamma{j}'*(Y{j}\gamma{j}); % Compute weighted norms
    phi(j)= phi_j;
    
    % compute BETA 
    beta(j)= log(det(Y{j})) + ngamma(j) - phi_j*log_C2;
end
    

%%
% 
% % possible number of outliers- array from 0 to the association with less outliers
% Noutliers_possible= unique(T(:,end));
% 
% % Generate the variables: h, H, gamma, ngamma
% for Noutliers_ind= 1:length(Noutliers_possible)
%     Noutliers= Noutliers_possible(Noutliers_ind);
%     if Noutliers == Nz, disp('all outliers'); return, end;
%     
%     ind= find(T(:,end) == Noutliers);
%     T_temporal= T(ind,:);
%     psi= length(ind);
%     
%     % make R a block diagonal with higher dimension
%     R= kron(eye(Nz-Noutliers),R_2x2);
% 
%     % stuck the models in different orders
%     H= cell(psi,1);
%     h= cell(psi,1);
%     gamma= cell(psi,1);
%     ngamma= zeros(psi,1);
%     Y= cell(psi,1);
%     for j= 1:psi
%         for i= 1:Nz
%             if T_temporal(j,i) ~= 0
%                 H{j}= [H{j};Hlm{T_temporal(j,i)}];
%                 h{j}= [h{j};hlm{T_temporal(j,i)}];
%             end
%         end
%         
%         % If there is an outlier eliminate the measurement associated with it
%         zj= z(:, T_temporal(j,1:Nz) ~= 0 );
%         
%         gamma{j}= zj(:) - h{j};
%         Y{j}= R + H{j}*P*H{j}';
%         
%         ngamma(j)= gamma{j}'*(Y{j}\gamma{j}); % Compute weighted norms      
%     end
%     
%     if any( ngamma < chi2inv(0.99, (Nz-Noutliers)*2) )
%         T= T_temporal;
%         break;
%     end;
% end
%     
    
%     
% % Select only the associations with less outliers
% for Noutliers= 0:Nz
%     ind= find(T(:,Nz+1) == Noutliers);
%     
%     if ~isempty(ind)
%         T= T(ind,:);
%         psi= length(ind);
%         
%         break
%     end
% end
% 
