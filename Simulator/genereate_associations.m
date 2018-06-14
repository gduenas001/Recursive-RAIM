
function genereate_associations(z,R_2x2,GATE)

global lm T hlm Hlm psi
global XX PX

Nxv= 3; % number of vehicle pose states
dz= 2;  % d.o.f. of each measurement
Nlm= size(lm,2); % number of features already in map
Nz= size(z,2); % number of measurements

% Generate individual models to compute NIS and gate later
hlm= cell(Nlm,1);
Hlm= cell(Nlm,1);
Slm= cell(Nlm,1);
for l= 1:Nlm
    [hlm{l},Hlm{l}]= observe_model_localization(XX,lm(:,l));
    Slm{l}= Hlm{l}*PX*Hlm{l}' + R_2x2;
end

% Compute the NIS distances for each possible individual association
nis= zeros(Nz,Nlm);
for i= 1:Nz
    for l= 1:Nlm
        v= z(:,i) - hlm{l};
        nis(i,l)= v'*(Slm{l}\v);
    end
end

% Eliminate very improbable associations using individual gates
A= cell(Nz,1);
for i= 1:Nz
    A{i}= [0, find(nis(i,:) < GATE)];
end

% Generate table with all possible combinations
T= allcomb(A{:}); psi= size(T,1); 
T= [T,-ones(psi,1)]; T(1,end)= Nz; % Add a column with the number of outliers

% Eliminate impossible associations, i.e. repeated landmarks
j= 2;
while j <= psi
    theta= sort(T(j,:));
    isoutlier= theta == 0;
    Noutliers= sum(isoutlier);
    if Noutliers
        ind= find(isoutlier,1,'last') + 1;
    else
        ind= 1;
    end
    
    isrep= ~all(diff(theta(ind:end)));    
    if isrep
        T(j,:)= [];
        psi= psi - 1;
    else
        T(j,Nz+1)= Noutliers;
        j= j + 1;
    end
end


