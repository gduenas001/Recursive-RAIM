function [Vn,Gn]= add_control_noise(G)
%function [V,G]= add_control_noise(V,G,Q, addnoise)
%
% Add random noise to nominal control values. We assume Q is diagonal.

global PARAMS SWITCH

if SWITCH.control_noise
    Vn= PARAMS.v + randn(1)*sqrt(PARAMS.Q(1,1));
    Gn= G + randn(1)*sqrt(PARAMS.Q(2,2));
end

