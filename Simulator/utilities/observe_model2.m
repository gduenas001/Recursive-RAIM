
function [h,Hz,Hx]= observe_model2(x,z,idf)

% fpos= idf*2 - 1; % position of xf in state

r= z(1,idf);
theta= z(2,idf);
c= cos(x(3)+theta);
s= sin(x(3)+theta);

h= [x(1) + r*c;
    x(2) + r*s];

Hz= [c, -r*s;
     s, r*c];
 
Hx= [1, 0, -r*s;
     0, 1,  r*c];




 
 







