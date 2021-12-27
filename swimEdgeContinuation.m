%SWIMEDGECONTINUATION Computes a fold of SwIM with particular
% by numerical continuation.
% Inputs:
% v0 - swimmer speed
% alpha - swimmer shape
% flow - fluid velocity field struct
% q_eq - swimming fixed point [x,y,theta]
% t1,t2 - tangent vectors to SwIM at q_eq
% beta0 - initial guess for initial condition parameter in SwIM
% Ti - initial integration time
% Tf - final integration time
% dT - step in integration time (optional)
% dobj_sign - desired sign of the objective function derivative (optional,
% default = -1)
%
% Outputs:
% betac - vector of beta coordinates on fold 
% Xfc - x coordinate on fold
% Yfc - y coordinate on fold
% THfc - theta coordinate on fold
% nfc - normal vector to SwIM along fold
% fvalc,exitFlagc, outputc, dobjc - fsolve outputs
% Tc - vector of integration time coorinates on fold
function [betac,Xfc,Yfc,THfc,nfc,fvalc,exitFlagc,outputc,dobjc,Tc] = swimEdgeContinuation(v0,alpha,flow,q_eq,t1,t2,beta0,Ti,Tf,dT,dobj_sign)
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
    if exitFlag > 0 && -dobj_sign*dobj < 0% record solution if fsolve succeeded and found a solution with negative slope (sign must be consistent with the slope near q_eq)
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