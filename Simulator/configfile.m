%%% Configuration file
%%% Permits various adjustments to parameters of the SLAM algorithm.
%%% See ekfslam_sim.m for more information

addpath('./utilities');
set(0,'DefaultFigureWindowStyle','normal')
format compact

% Initialise states and other global variables
global XX PX BIAS GAMMA_D LM POSES  DATA PARAMS SWITCH IA MSMTS LB LM_obstacle

% % Load table for the lower bounds of the non-centrality parameter
% load('poses.mat'); POSES= poses; clear poses;
% load('LM.mat'); LM= LM'; 


% Tunel LMs
POSES= [0,0,0; 10000,0,0];
posLM(1,:)= linspace(0,3000,102);
negLM(1,:)= linspace(0,3000,102);

% well spaced
posLM(2,1:7)= linspace(-15,-15,7);
negLM(2,1:7)= linspace(15,15,7);

% Bad spaced
posLM(2,8:102)= linspace(-2.5,-2.5,102-7);
negLM(2,8:102)= linspace(2.5,2.5,102-7);
LM= [posLM,negLM];




% Number of time epochs to run the simulation
PARAMS.numSteps= 1300;
PARAMS.dt= 0.1;

% Initial control inputs
iwp= 1;
G= deg2rad(0);

% control parameters
PARAMS.m= 3;
PARAMS.at_waypoint= 5;
PARAMS.V= 10;
PARAMS.maxRateG= deg2rad(50);
PARAMS.maxG= deg2rad(30);
PARAMS.wheelbase= 1; % metres, vehicle wheel-base
PARAMS.veh= [0 -4    -4   0; 
              0 0.5  -0.5 0]; % vehicle animation

% control noises
PARAMS.sigmaV= 0.1; % m/s ( default= 0.3 )
PARAMS.sigmaG= deg2rad(1); % radians ( default= (3.0*pi/180) )
PARAMS.Q= [PARAMS.sigmaV^2, 0; 0, PARAMS.sigmaG^2];

% observation parameters
PARAMS.maxRange= 25;
PARAMS.m_F= 2; % d.o.f. of one measurement

% observation noises
PARAMS.sigmaR= 0.3; % metres ( default 0.1 )
PARAMS.sigmaB= deg2rad(2); % radians ( default (1.0*pi/180) )
PARAMS.R= [PARAMS.sigmaR^2 0; 0 PARAMS.sigmaB^2];

% Integrity 
PARAMS.I_y= 1e-10;
PARAMS.I_T= 0.01;
PARAMS.I_FOV= 1e-9;
PARAMS.I_REQ= 1e-5;
PARAMS.alert_limit= 1;
PARAMS.T2= chi2inv(1-PARAMS.I_T,PARAMS.m_F);
PARAMS.P_IA_max= PARAMS.I_FOV*4;
PARAMS.n_W= 10; % Window of epochs
PARAMS.C_REQ= 1e-7;
PARAMS.Preceding_Horizon= 3; % epochs
PARAMS.Epoch= 0; 
PARAMS.P_UA= 10^-3; % assuming that it is constant for the whole landmarks.
PARAMS.P_ME= 0; % Misdetection probability

% switches
SWITCH.control_noise= 1; % if 0, velocity and gamma are perfect
SWITCH.sensor_noise= 1; % if 0, measurements are perfect
SWITCH.inflate_noise= 0; % if 1, the estimated Q and R are inflated (ie, add stabilising noise)
SWITCH.heading_known= 0; % if 1, the vehicle heading is observed directly at each iteration
SWITCH.batch_update= 1; % if 1, process scan in batch, if 0, process sequentially
SWITCH.seed_random= 3; % if not 0, seed the randn() with its value at beginning of simulation (for repeatability)
SWITCH.use_IEKF= 0; % if 1, use iterated EKF for updates, if 0, use normal EKF
SWITCH.profile= 0; % if 1, turn on MatLab profiling to measure time consumed by simulator functions
SWITCH.graphics= 1; % if 0, avoids plotting most animation data to maximise simulation speed
SWITCH.update_global= 0; % if 1, This alternative is a "global constraint" model devised by Jose Guivant, and may have better linearisation properties than the conventional range-bearing model.
SWITCH.association= 1; % if 0, associations are given; if 1, they are estimated using the LNN
SWITCH.consistency= 0; % If 1, checks innovation vector consistency.
SWITCH.ME= 0; % If 0, no mis-extractions; if 1, there are mis-extractions.
SWITCH.update_if_necessary= 0; % If 1, only updates when integrity improves.
if SWITCH.seed_random, rand('state',SWITCH.seed_random), randn('state',SWITCH.seed_random), end
if SWITCH.profile, profile on -detail builtin, end

%% INITIALIZATIONS 
 
% True & estimated state
xtrue= POSES(1,:)';
XX= xtrue; 
BIAS= zeros(length(xtrue),1);
GAMMA_D= 0;

% initial pose covariance
PX= [eps,   0,   0;
       0, eps,   0;
       0,   0, eps];


PARAMS.ftag= 1:size(LM,2);     % identifier for each landmark
da_table= zeros(1,size(LM,2)); % data association table 

% more initializations
step= 0; IA= 0;
DATA.epsXX= zeros(5000,3);
DATA.stdXX= zeros(5000,3);
DATA.eps= zeros(5000,1);
DATA.stdEps= zeros(5000,1);
DATA.P_HMI= zeros(5000,1);
DATA.P_HMI_CA= zeros(5000,1);
DATA.path= zeros(5000,2);
DATA.PnMA_k= ones(5000,1);
DATA.PnMA_K= ones(5000,1);
DATA.IA= zeros(5000,1);
DATA.numAssoc= zeros(5000,1);

DATA.gamma= cell(5000,1);
DATA.stngamma= cell(5000,1);
DATA.PCAt= ones(5000,1);
DATA.realPCA= zeros(5000,1);
DATA.calcPCA= zeros(5000,1);
DATA.bias= zeros(5000,3);
DATA.bias_interest= zeros(5000,1);
DATA.T_RB= zeros(5000,1);
DATA.q_RB= zeros(5000,1);
DATA.lambda2= zeros(5000,1);
DATA.lambda2_current= zeros(5000,1);



