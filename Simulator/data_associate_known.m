function [idft]= data_associate_known(z, idft, EFOV)
%function [zf,idf,zn, table]= data_associate_known(x,z,idz, table)

global XX LM hlm Hlm PARAMS

if PARAMS.P_ME ~= 0, error('Probability of ME has to be 0 for the known DA'), end;

Nz= size(z,2);
n_L= size(LM,2);

% Initialize lm models to be filled in the next function
hlm= cell(n_L,1);
Hlm= cell(n_L,1);

% Get all visible landmarks, assuming no mis-extractions here
lm_ind= get_visible_landmarks(XX,PARAMS.maxRange+EFOV, 0);


for i= 1:Nz
    for l= 1:n_L
        
        % It's in the extended FV and not removed by LS
        if ismember(l,lm_ind) 
            lm_model(l);
        end
        
    end
end




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function lm_model(lm_id)
% return normalised innovation squared (ie, Mahalanobis distance) and normalised distance

global XX LM hlm Hlm



% auxiliary values
dx= LM(1,lm_id) - XX(1); 
dy= LM(2,lm_id) - XX(2);
d2= dx^2 + dy^2;
d= sqrt(d2);

xd= dx/d;
yd= dy/d;
xd2= dx/d2;
yd2= dy/d2;

% predict z
h= [d;
    atan2(dy,dx) - XX(3)];

% calculate H
H = [-xd -yd 0; 
      yd2 -xd2 -1];

% Store the values in the global varibles
hlm{lm_id}= h;
Hlm{lm_id}= H;










