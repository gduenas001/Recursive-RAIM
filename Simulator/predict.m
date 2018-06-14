function [xtrue,xx,Px,Gx]= predict (xtrue, xx, Px, G)

global BIAS PARAMS 

% True State
xtrue= [xtrue(1) + PARAMS.V*PARAMS.dt*cos(G+xtrue(3));
        xtrue(2) + PARAMS.V*PARAMS.dt*sin(G+xtrue(3));
        pi_to_pi( xtrue(3)+ PARAMS.V*PARAMS.dt*sin(G) / PARAMS.wheelbase )];


    
% Add noise to the controls
Vn= PARAMS.V + normrnd(0, PARAMS.sigmaV); 
Gn= G + normrnd(0, PARAMS.sigmaG); 
% Vn= PARAMS.V;
% Gn= G;

% compute variables
s= sin(Gn + xx(3)); 
c= cos(Gn + xx(3));
vts= Vn*PARAMS.dt*s; 
vtc= Vn*PARAMS.dt*c;

% Predict state
xx= [xx(1) + vtc;
     xx(2) + vts;
     pi_to_pi( xx(3)+ Vn*PARAMS.dt*sin(Gn) / PARAMS.wheelbase )];

% jacobians   
Gx= [1 0 -vts;
     0 1  vtc;
     0 0 1];
Gu= [PARAMS.dt*c,                       -vts;
     PARAMS.dt*s,                        vtc;
     PARAMS.dt*sin(Gn)/PARAMS.wheelbase, Vn*PARAMS.dt*cos(Gn)/PARAMS.wheelbase];

% predict covariance
Px= Gx*Px*Gx' + Gu*PARAMS.Q*Gu';


% Bias prediction
BIAS= Gx * BIAS;
