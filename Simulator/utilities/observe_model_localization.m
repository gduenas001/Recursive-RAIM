function [z,H]= observe_model_localization(x, lmi)
%function [z,H]= observe_model(x, idf)
%
% INPUTS:
%   x - state vector
%   idf - index of feature order in state
%
% OUTPUTS:
%   z - predicted observation
%   H - observation Jacobian
%
% Given a feature index (ie, the order of the feature in the state vector),
% predict the expected range-bearing observation of this feature and its Jacobian.
%
% Tim Bailey 2004.

% auxiliary values
dx= lmi(1)  -x(1); 
dy= lmi(2)-x(2);
d2= dx^2 + dy^2;
d= sqrt(d2);
xd= dx/d;
yd= dy/d;
xd2= dx/d2;
yd2= dy/d2;

% predict z
z= [d;
    atan2(dy,dx) - x(3)];

% calculate H
H(:,1:3)        = [-xd -yd 0; 
                    yd2 -xd2 -1];
