
function [idf]= data_associate_LNN (z)

global hlm Hlm PARAMS LM

Nz= size(z,2);
Nlm= size(LM,2);

% Initialize associations
idf= zeros(1,Nz);

% Initialize lm models to be filled in the next function
hlm= cell(Nlm,1);
Hlm= cell(Nlm,1);

% Create T associations table where the each columns corresponds to a lm
IIVN= associations_tables(z, Nz, Nlm, PARAMS.R, PARAMS.gate);


% Choose best validated association for each feature
[minIIN_value, minIIN_ind]= min(IIVN,[],2);
ind_validated= minIIN_value < PARAMS.gate;
idf(ind_validated)= minIIN_ind(ind_validated);




% % Best association by rows (by feature)
% for i=1:Nz
%     validated_assoc= IIVN(i,:) ~= 0;
%     if sum(validated_assoc) > 1
%         [chosenValue,chosenInd]= min(IIVN(i,validated_assoc));
%         inds= find(validated_assoc);
%         chosenInd= inds(chosenInd);
%         IIVN(i,validated_assoc)= 0;
%         IIVN(i,chosenInd)= chosenValue;
%     end
% end
% 
% % Best association by columns (by lm)
% for l=1:Nlm
%     validated_assoc= IIVN(:,l) ~= 0;
%     sum_validated_assoc= sum(validated_assoc);
%     if sum_validated_assoc == 1
%         idf(find(validated_assoc))= l;
%         if sum_validated_assoc > 1
%             [chosenValue,chosenInd]= min(IIVN(validated_assoc,l));
%             inds= find(validated_assoc);
%             chosenInd= inds(chosenInd);
%             idf(chosenInd)= l;
%             
%             IIVN(validated_assoc,l)= 0;
%             IIVN(chosenInd,l)= chosenValue;
%         end
%     end
% end
% 




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [IIVN]= associations_tables(z, Nz, Nlm, R, gate)

% Create the nis table
IIVN= ones(Nz,Nlm)*gate;
for i= 1:Nz
    for l= 1:Nlm
        [IIVN(i,l)]= lm_model_and_nis(z(:,i),R,l);        
    end
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
function [nis]= lm_model_and_nis(z,R,lm_id)
% return normalised innovation squared (ie, Mahalanobis distance) and normalised distance

global XX PX LM hlm Hlm


% auxiliary values
dx= LM(1,lm_id)  -XX(1); 
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


% Innovation vector
v= z-h; v(2)= pi_to_pi(v(2));

% Quick check to discard associations
if abs(v(1)) > 10 || abs(v(2)) > deg2rad(45), nis= inf; return, end;

% Compute IIVN 
S= H*PX*H' + R; 

nis= v'/S*v;
































