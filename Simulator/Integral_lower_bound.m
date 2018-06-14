

%% Vary degrees of freedom: p
clear; clc; close all;

p= 2:8:60;
lambda2= 4; lambda= sqrt(lambda2);

for i= 1:length(p)
    gamma02= chi2rnd(p(i),[100000,1]);
    gamma0= sqrt(gamma02);

    fun= gamma0 - ((p(i)-1)/2) * (log(gamma0) - log(lambda))./(gamma0 - lambda);
    
    subplot(2,round(length(p)/2),i); hold on; grid on;
    title(['lambda^2: ',num2str(lambda2), '    d.o.f.: ', num2str(p(i))]); 
    histogram(fun,'Normalization','pdf');
    axis([-2,6,0,0.6]);
end


%% Vary noncetrality parameter: lambda2
clear; clc; close all;

lambda2= 2:8:48; lambda= sqrt(lambda2);
p= 6;

for j= 1:length(lambda)
    gamma02= chi2rnd(p,[100000,1]);
    gamma0= sqrt(gamma02);

%     fun= gamma0 - ((p-1)/2) * (log(gamma0) - log(lambda(j)))./(gamma0 - lambda(j));
    fun2=  (log(gamma0) - log(lambda(j)))./(gamma0 - lambda(j));
%     (p-1)/2
    subplot(2,round(length(lambda)/2),j); hold on; grid on;
    title(['lambda^2: ',num2str(lambda2(j)), '    d.o.f.: ', num2str(p)]); 
    histogram(fun2,'Normalization','pdf');
%     histogram(gamma0,'Normalization','pdf');
%     axis([-2,6,0,0.6]);
end
















