function [G,iwp]= compute_steering(xx, iwp, G)
%function [G,iwp]= compute_steering(x, wp, iwp, minD, G, rateG, maxG, dt)
%
% INPUTS:
%   x - true position
%   wp - waypoints
%   iwp - index to current waypoint
%   minD - minimum distance to current waypoint before switching to next
%   G - current steering angle
%   rateG - max steering rate (rad/s)
%   maxG - max steering angle (rad)
%   dt - timestep
%
% OUTPUTS:
%   G - new current steering angle
%   iwp - new current waypoint

global PARAMS POSES


% determine if current waypoint reached
while 1
    cwp= POSES(iwp,:);
    d2= (cwp(1)-xx(1))^2 + (cwp(2)-xx(2))^2;
    
    if d2 < PARAMS.at_waypoint^2 
        iwp= iwp+1; % switch to next
        if iwp > size(POSES,1) % reached final waypoint, flag and return
            iwp=0;
            return;
        end
        cwp= POSES(iwp,:); % next waypoint
    else
        break;
    end
end

% compute change in G to point towards current waypoint
deltaG= pi_to_pi(atan2(cwp(2)-xx(2), cwp(1)-xx(1)) - xx(3) - G);

% limit rate
maxDelta= PARAMS.maxRateG*PARAMS.dt;
if abs(deltaG) > maxDelta
    deltaG= sign(deltaG)*maxDelta;
end

% limit angle
G= G+deltaG;
if abs(G) > PARAMS.maxG
    G= sign(G)*PARAMS.maxG;
end




