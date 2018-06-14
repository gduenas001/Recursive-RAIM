
function do_plots(xtrue, z, h, step)

global XX PARAMS SWITCH DATA
persistent pcount

% Plots
if SWITCH.graphics
    
    % True vehicle
    xt= transformtoglobal(PARAMS.veh, xtrue);
    set(h.xt, 'xdata', xt(1,:), 'ydata', xt(2,:))
    
    % Estimate
    xv= transformtoglobal(PARAMS.veh, XX);
    set(h.xv, 'xdata', xv(1,:), 'ydata', xv(2,:))
        
    % Covariance estimate
    pvcov= make_vehicle_covariance_ellipse();
    set(h.vcov, 'xdata', pvcov(1,:), 'ydata', pvcov(2,:))
    
    % Plot the path
    if isempty(pcount), pcount= 1; else pcount= pcount+1; end;
    
%     if pcount == 5 % infrequently
%         pcount=0;
%         set(h.path_true, 'xdata', POSES(1:step,1), 'ydata', POSES(1:step,2));
%         set(h.path, 'xdata', DATA.path(2:step,1), 'ydata', DATA.path(2:step,2));
%     end
    
    % plots related to observations
    if ~isempty(z) 
        plines= make_laser_lines (z,xtrue);
        set(h.obs, 'xdata', plines(1,:), 'ydata', plines(2,:))
    end
    
    % Plot the P(HMI) online
    if DATA.P_HMI(step) > h.figHMI.CurrentAxes.YLim(2)
        h.figHMI.CurrentAxes.YLim(2)= DATA.P_HMI(step) + 0.2*DATA.P_HMI(step);
    end
    set(h.HMI, 'xdata', [1:step], 'ydata', DATA.P_HMI(1:step));
    set(h.HMI_CA, 'xdata', [1:step], 'ydata', DATA.P_HMI_CA(1:step));
        
    % Plot the P(MA)_k online (at current time)
    set(h.PMA_k, 'xdata', [1:step], 'ydata', 1 - DATA.PnMA_k(1:step));
    set(h.PMA_K, 'xdata', [1:step], 'ydata', 1 - DATA.PnMA_K(1:step));
    
    % Plot epsilon error online
    if DATA.stdEps(step) > h.figEps.CurrentAxes.YLim(2)
        h.figEps.CurrentAxes.YLim(2)= DATA.stdEps(step) + 0.2*DATA.stdEps(step);
    end
    set(h.eps, 'xdata', [1:step], 'ydata', abs(DATA.eps(1:step)));
    set(h.cov, 'xdata', [1:step], 'ydata', DATA.stdEps(1:step));    
    
    % Plot the IAs online
    set(h.IA, 'xdata', [1:step], 'ydata', DATA.IA(1:step));
    set(h.numAssoc, 'xdata', [1:step], 'ydata', DATA.numAssoc(1:step)); 

    % Plot the bias online
    set(h.BIAS, 'xdata', [1:step], 'ydata', DATA.bias_interest(1:step));

    % Plot the detector online
    set(h.q_RB, 'xdata', [1:step], 'ydata', DATA.q_RB(1:step));
    set(h.T_RB, 'xdata', [1:step], 'ydata', DATA.T_RB(1:step));

    
    drawnow
end










