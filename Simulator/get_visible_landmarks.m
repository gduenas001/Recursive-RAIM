
function [idf]= get_visible_landmarks(xtrue,maxRange,isMD)

global PARAMS LM

% Select set of landmarks that are visible within vehicle's semi-circular field-of-view
dx= LM(1,:) - xtrue(1);
dy= LM(2,:) - xtrue(2);

% incremental tests for bounding semi-circle
ii= find(abs(dx) < maxRange & abs(dy) < maxRange ... % bounding box
      & (dx.^2 + dy.^2) < maxRange^2);           % bounding circlew
  %       & (dx*cos(phi) + dy*sin(phi)) > 0 ...  % bounding line

  
% Eliminate the non-detected lm's using the probabitlity of detection
if isMD
    indMD= boolean( binornd(1, PARAMS.P_ME, [1,length(ii)]) );
    ii(indMD)= [];
end

idf= PARAMS.ftag(ii);

