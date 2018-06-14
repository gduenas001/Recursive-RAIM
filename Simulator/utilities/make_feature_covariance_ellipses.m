function p= make_feature_covariance_ellipses(x,P)
% compute ellipses for plotting feature covariances
N= 10;
inc= 2*pi/N;
phi= 0:inc:2*pi;
circ= 2*[cos(phi); sin(phi)];

lenx= length(x);
lenf= (lenx-3)/2;
p= zeros (2, lenf*(N+2));

ctr= 1;
for i=1:lenf
    ii= ctr:(ctr+N+1);
    jj= 2+2*i; jj= jj:jj+1;
    
    p(:,ii)= make_ellipse(x(jj), P(jj,jj), circ);
    ctr= ctr+N+2;
end

%
%

