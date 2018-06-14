
function [PHMI_CA,Px]= PHMI_with_CA(xx, Px, l)

% if issymmetric(P)
%     PHMI_CA= 1 - mvncdf([l; l], [0;0], P);
% else
%     P= make_sym(P);
%     PHMI_CA= 1 - mvncdf([l; l], [0;0], P);
% end
% alpha= [cos(xx(3)) ; sin(xx(3))];
alpha= [-sin(xx(3)) ; cos(xx(3))];
PHMI_CA= 2*normcdf(-l,0, sqrt(alpha' * Px * alpha));

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% function M= make_sym(M)
% M= (M + M')/2;




