function p= make_vehicle_covariance_ellipse()

global XX PX

% compute ellipses for plotting vehicle covariances
N= 10;
inc= 2*pi/N;
phi= 0:inc:2*pi;
circ= 3*[cos(phi); sin(phi)];

p= make_ellipse(XX(1:2), PX(1:2,1:2), circ);

