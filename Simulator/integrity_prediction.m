
function [xPrev,ensureTime]= integrity_prediction(xPrev, Px, G, iwp, PnMA0)

global PARAMS 

ensureSteps= 0;
while 1
    [G, iwp]= compute_steering(xPrev, iwp, G);
    [~, xx, Px]= predict(xPrev,xPrev,Px,G);
    
    PHMI_nMA= PHMI_with_CA(xx,Px(1:2,1:2),PARAMS.alert_limit);
    PHMI= 1 + ( PHMI_nMA - 1 )* PnMA0 ;
    
    if PHMI > PARAMS.I_REQ
        ensureTime= ensureSteps * PARAMS.dt;
        break;
    end
    
    xPrev= xx;
    ensureSteps= ensureSteps + 1;
end

% PnMA= prob_nMA(xx,Px);









































