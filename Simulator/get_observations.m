function [z,idf]= get_observations(xtrue)
% INPUTS:
%   x - vehicle pose [x;y;phi]
%   lm - set of all landmarks
%   idf - index tags for each landmark
%   rmax - maximum range of range-bearing sensor 
%
% OUTPUTS:
%   z - set of range-bearing observations
%   idf - landmark index tag for each observation

global LM PARAMS SWITCH

idf= get_visible_landmarks(xtrue,PARAMS.maxRange,SWITCH.ME);

% Select the visible landmarks
lm= LM(:,idf);

% Compute exact observation
dx= lm(1,:) - xtrue(1);
dy= lm(2,:) - xtrue(2);

z= [sqrt(dx.^2 + dy.^2);
    atan2(dy,dx) - xtrue(3)];
    

