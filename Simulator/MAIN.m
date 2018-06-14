
dbstop if error
dbclear if error


clear; close all;
configfile; 

h= setup_animations();


% *****************    MAIN LOOP    *****************
for step= 2:PARAMS.numSteps
    disp(['Step: ',num2str(step)]);
    
    % Compute controls
    [G,iwp]= compute_steering(xtrue, iwp, G);
    
    % EKF predict step
    [xtrue,XX,PX,Gx]= predict (xtrue,XX,PX,G);
    
    % Integrity Monitoring    
    [P_HMI_worst,F_mag,Fault_slope,H,L,L_p,L_pp,Y]= integrity_monitoring_fault...
        (Px,G,H,L,L_p,L_pp,A_k,Y,T,P_CA);
    
    % Get measurements
    [z,idft]= get_observations(xtrue);
    z= add_observation_noise(z);
    
    % DA
    if  ~isempty(z)
        [idf,DATA.numAssoc(step)]= data_associate_LNN_LS(z, remove_lm_ind, EFV);
                
        % Store associations data
        DATA.IA(step)= any( (idft - idf).*idf );
    else
        idf= []; 
        DATA.numAssoc(step)= 0;
        DATA.IA(step)= 0;
    end
    
    [q_RB, T_RB]= EKF_update(z, idf, step);
    
    % Store DATA
    store_data(step, P_HMI_nMA_UA, PnMA_k, xtrue, q_RB, T_RB);
        
    % Plots    
    do_plots(xtrue, z, h, step);

%     pause(0.1);
end 
% *****************  END OF MAIN LOOP    *****************

% post_processing_and_plots(step)



