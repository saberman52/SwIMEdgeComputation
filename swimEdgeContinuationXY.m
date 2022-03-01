%SWIMEDGECONTINUATIONXY Computes a frontier---fold of SwIM with particular
%orientation---by numerical continuation.
%Uses spacing between SwIM edge points as a criterion for keeping a step of
%a certain size or reducing the step size.
% Inputs:
% v0 - swimmer speed
% alpha - swimmer shape
% flow - fluid velocity
% q_eq - swimming fixed point [x,y,theta]
% t1,t2 - tangent vectors to SwIM at q_eq
% beta0 - initial guess for swimmer initial condition parameter
% Ti - initial integration time
% Tf - final integration time
% dT - step in integration time (optional)
% dobj_sign - desired sign of the objective function derivative (optional,
% default = -1)
%
% Outputs:
% betac - initial condition parameters leading to frontier
% Xfc - x coordinate on frontier
% Yfc - y coordinate on frontier
% THfc - theta coordinate on frontier
% nfc - normal vector to SwIM along frontier
% fvalc,exitFlagc, outputc, dobjc 0- fsolve outputs
% Tc - time coorinate on frontier
function [betac,Xfc,Yfc,THfc,nfc,fvalc,exitFlagc,outputc,dobjc,Tc] = swimEdgeContinuationXY(v0,alpha,flow,q_eq,t1,t2,beta0,Ti,Tf,dT,dobj_sign)
ep = 1e-3; % initial step size away from equilibrium
tols = 1e-9*[1,1]; % integration tolerances

if nargin < 10
    dT = 0.1; % initial change in integration time
end
if nargin < 11
    dobj_sign = -1;
end
% vectors to store along continuation
Tc = []; 
betac = [];
Xfc = [];
Yfc = [];
THfc = [];
nfc = [];
fvalc = [];
exitFlagc = [];
outputc = [];
dobjc = [];

T = Ti; % current integration time

while T < Tf
    % try to find Swim edge point with these parameters
    [beta,Xf,Yf,THf,nf,fval,exitFlag,output,dobj] = swimEdge0(beta0,T,ep,q_eq,t1,t2,v0,alpha,flow,tols);
    disp(['exitFlag = ' num2str(exitFlag)])
    % check spacing between subsequent points in xy plane
    if T > Ti
        close_enough = sqrt((Xf - Xfc(end)).^2 + (Yf - Yfc(end)).^2) < 0.05;
    else
        close_enough = true;
    end
    if exitFlag > 0 && -dobj_sign*dobj < 0 && close_enough% if fsolve succeeded and found a solution with negative slope (i.e. a frontier, not a positive-slope fold)
           
            betac = [betac; beta];
            Xfc = [Xfc; Xf];
            Yfc = [Yfc; Yf];
            THfc = [THfc; THf];
            nfc = [nfc; nf];
            fvalc = [fvalc; fval];
            exitFlagc = [exitFlagc; exitFlag];
            outputc = [outputc; output];
            dobjc = [dobjc; dobj];
            Tc = [Tc; T];
    elseif dT > 1e-5 && T ~= Ti % arbitrary threshold for step size
        dT = dT/2;
        T = Told + dT; % retry with smaller step
        continue
    elseif dT > 1e-5 && T == Ti
        disp('Failed on initial guess.')
        disp(['exitFlag = ' num2str(exitFlag)])
        disp(['dobj = ' num2str(dobj)])
        break
    else
        disp('Terminated because dT fell below threshold.')
        break
    end
    % store initial parameters
    Told = T;
    % set new parameters
    T = T+dT; 
    beta0 = beta;
end

end