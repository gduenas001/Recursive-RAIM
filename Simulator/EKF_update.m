function [q_RB, T_RB]= EKF_update(z, idf, step)

global XX PX PARAMS hlm Hlm DATA

% Define some variables
n_L= DATA.numAssoc(step);
n= n_L * PARAMS.m_F;

% If no lms are associated --> return!
if n == 0
    q_RB= 0;
    T_RB= 0;
    return;
end

% Update detector threshold
T_RB= chi2inv(1 - PARAMS.C_REQ, n);

% Remove non-associated msmts
z(:, idf == 0)= [];
idf( idf== 0) = [];

% make R a block diagonal with higher dimension
R= kron(eye(n_L), PARAMS.R);

% create the models for the association
h= zeros(PARAMS.m_F * n_L,1);
H= zeros(PARAMS.m_F * n_L,3);
for i= 1:n_L
    ind= i*PARAMS.m_F - 1;
    h(ind:ind+1)= hlm{idf(i)};
    H(ind:ind+1,:)= Hlm{idf(i)};
end

% Compute innovations
gamma= z(:) - h;
gamma(2:2:end)= pi_to_pi(gamma(2:2:end));
Y= H*PX*H' + R;

% Save previous estimate
XX_bar= XX;
PX_bar= PX;

% Update the estimate
K= PX*H'/Y;
PX= PX - K*H*PX;
XX= XX + K*gamma;


%% DETECTOR & BIAS %%

% Create model
diffz= z(:) - h;
diffz(2:2:end)= pi_to_pi( diffz(2:2:end) );
z_star= diffz + H*XX_bar;
z_star(2:2:end)= pi_to_pi( z_star(2:2:end) );
l= [z_star; XX_bar];
D= [H; eye(3)];
r= l - D*XX;
r(2:2:n)= pi_to_pi( r(2:2:n) );
r(end)= pi_to_pi(r(end));

% Detector
q_RB= r(1:n)' / R * r(1:n) + r(n+1:end)' / PX_bar * r(n+1:end);

