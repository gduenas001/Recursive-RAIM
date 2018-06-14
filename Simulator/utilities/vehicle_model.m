function xv= vehicle_model(xv, G)
%
% INPUTS:
%   xv - vehicle pose [x;y;phi]
%   V - velocity
%   G - steer angle (gamma)
%   WB - wheelbase
%   dt - change in time
%
% OUTPUTS:
%   xv - new vehicle pose

global PARAMS

xv= [xv(1) + PARAMS.v*PARAMS.dt*cos(G+xv(3,:)); 
     xv(2) + PARAMS.v*PARAMS.dt*sin(G+xv(3,:));
     pi_to_pi(xv(3) + PARAMS.v*PARAMS.dt*sin(G)/PARAMS.wheelbase)];
 