%SWIMEDGECONTINUATIONCIRCDET Computes SwIM folds by numerical continuation in
%the (T,beta) plane, where T is the hitting time to the fold and beta is
%the initial condition parameter.
%Uses a zero finder with det dr/d(T,beta) as objective function.
%Initially search for point on the fold, and then continue from there
% Inputs:
% v0 - swimmer speed
% alpha - swimmer shape
% flow - fluid velocity
% q_eq - swimming fixed point [x,y,theta]
% t1,t2 - tangent vectors to SwIM at q_eq
% beta0 - initial guess for swimmer initial condition parameter
% T0 - initial integration time
% arcmax - maximum arclength of curve in (T,beta) plane
% dT - step size for time coordinate for second solution (optional)
% ep0 - initial step size in (T,beta) plane (optional)
% epMin - threshold for minimum step size (optional)
%
% Outputs:
% gammac - angle parameter in (T,beta) plane
% epc - step size
% betac - initial condition parameters leading to fold
% Xfc - x coordinate on fold
% Yfc - y coordinate on fold
% THfc - theta coordinate on fold
% nfc - normal vector to SwIM along fold
% fvalc,exitFlagc, outputc, dobjc 0- fsolve outputs
% Tc - time coorinate on fold
function [gammac,betac,Xfc,Yfc,THfc,nfc,fvalc,exitFlagc,outputc,dobjc,Tc] = swimEdgeContinuationCircDet(v0,alpha,flow,q_eq,t1,t2,beta0,T0,arcmax,dT,ep0,epMin)
ep = 1e-3; % initial step size away from equilibrium
tols = 1e-9*[1,1]; % integration tolerances

if nargin < 10
    dT = 0.05; % initial step size for time step
end

if nargin < 11
    ep0 = 0.1; % initial perturbation size
end

if nargin < 12
    epMin = 1e-8; % minimum threshold for perturbation step
end

gamma0 = 0; % initial guess for angle parameter
% vectors to store along continuation
gammac = [];
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

arc = 0; % arclength of curve in (T,beta) plane

%% first two iterations: find initial point on fold and figure out gamma
for i = 1:2
    [beta,Xf,Yf,THf,nf,fval,exitFlag,output,dobj] = swimEdge0_det_fixedT(beta0,T0,ep,q_eq,t1,t2,v0,alpha,flow,tols);
    if exitFlag <= 0
        if i == 1
            disp('Failed on initial guess.')
        else
            disp('Failed on second iteraion')
        end
        disp(['exitFlag = ' num2str(exitFlag)])
        disp(['dobj = ' num2str(dobj)])
        return
    else
        betac = [betac; beta];
        Xfc = [Xfc; Xf];
        Yfc = [Yfc; Yf];
        THfc = [THfc; THf];
        nfc = [nfc; nf];
        fvalc = [fvalc; fval];
        exitFlagc = [exitFlagc; exitFlag];
        outputc = [outputc; output];
        dobjc = [dobjc; dobj];
        Tc = [Tc; T0];
        % set new parameters
        if i == 1
            T0 = Tc(end) + dT;
        else
            T0 = Tc(end);
        end
        beta0 = betac(end);
    end
end
% set initial gamma
gamma0 = atan2(betac(2)-betac(1),Tc(2)-Tc(1));
%% loop
while arc < arcmax
    % try to find Swim edge point with these parameters
    [gamma,T,beta,Xf,Yf,THf,nf,fval,exitFlag,output,dobj] = swimEdge0_det(gamma0,T0,beta0,ep,ep0,q_eq,t1,t2,v0,alpha,flow,tols);
    disp(['exitFlag = ' num2str(exitFlag)])
    if arc == 0
        gammaCompare = gamma0;
    else
        gammaCompare = gammac(end);
    end
    if exitFlag > 0 %&& abs(gamma - gammaCompare) < 0.1 % if fsolve succeeded and found a solution
        gammac = [gammac; gamma];
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
        % update arc length
        arc = arc + ep0;
        disp(['arc = ' num2str(arc) ', ep0 = ' num2str(ep0)])
    elseif ep0/2 > epMin %threshold for step size
        ep0 = ep0/2; % retry with smaller step
        disp(['Zero search failed, retrying with ep0 = ' num2str(ep0)])
        continue
    else
        disp('Terminated because ep0 fell below threshold.')
        break
    end
    % increase step up to a threshold if 
    % a) converged to previous solution sufficiently fast &&
    % b) ep0 is not too large &&
    % c) the difference between current and previous gamma is sufficiently
    % small
    if output.iterations < 10 && 1.2*ep0 < 0.2 && abs(gamma - gammac(end)) < 1e-3
        ep0 = 1.2*ep0;
    end
    % set new parameters
    gamma0 = gammac(end);
    T0 = Tc(end);
    beta0 = betac(end);
end

end