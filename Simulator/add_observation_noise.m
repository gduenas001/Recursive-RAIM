function z= add_observation_noise(z)
%function z= add_observation_noise(z,R, addnoise)
%
% Add random measurement noise. We assume R is diagonal.

global PARAMS SWITCH

if SWITCH.sensor_noise
    len= size(z,2);
    if len > 0
        z(1,:)= z(1,:) + randn(1,len)*sqrt(PARAMS.R(1,1));
        z(2,:)= z(2,:) + randn(1,len)*sqrt(PARAMS.R(2,2));
    end
end
