% SWIMEDGE0_DET searches for the initial condition and hitting time of a 
% point hitting the SwIM edge at time T, near a known point on the SwIM
% edge.
% Initial condition of form q0 = q_eq + ep*(cos(beta)*t1 +
% sin(beta)*t2), where (T,beta) = (T0,beta0) + ep0*(cos(gamma),sin(gamma)).
% Objective function = det(d r/ d(T,\beta))
% Inputs:
% gamma0 - initial guess for time/intial condition parameter
% T0 - known hitting time for the initial condition beta0
% beta0 - known initial condition parameter
% ep - initial step size away from equilibrium
% ep0 - initial step size away from known point
% q_eq - equilibrium point (row vector)
% t1 - one of the tangent vectors (row vector)
% t2 - other tangent vector
% v0 - swim speed
% alpha - swimmer shape parameter
% flow - fluid velocity field
% tols - numerical integration tolerances
%
% Outputs:
% gamma - final initial condition parameter
% T - final hitting time
% beta - final inital condition
% Xf - final x coordinate on SwIM edge
% Yf - final y coordinate
% THf - final theta coordinate
% nf - final normal vector to manifold (normalized)
% dobj - derivative of the objective function at the solution
function [gamma,T,beta,Xf,Yf,THf,nf,fval,exitFlag,output,dobj] = swimEdge0_det(gamma0,T0,beta0,ep,ep0,q_eq,t1,t2,v0,alpha,flow,tols)
%% compute intial condition of trajectory going to SwIM edge in time T
%% fsolve
% opt = optimoptions('fsolve','Display','iter-detailed','OptimalityTolerance',1e-8);
% [gamma,fval,exitFlag,output,dobj] = fsolve(@objective,gamma0,opt);
%% fzero
opt = optimset('Display','iter');
[gamma,fval,exitFlag,output] = fzero(@objective,gamma0,opt);
dobj = 0;

T = T0 + ep0*cos(gamma);
beta = beta0 + ep0*sin(gamma);
%% compute the resulting trajectory
q00 = q_eq + ep*(cos(beta)*t1 + sin(beta)*t2); % initial condition
% evaluate trajectory to get final points on SwIM edge
[X,Y,TH,J,~,~] = swimmerTangentU_Jt(q00(1),q00(2),q00(3),[0 T],flow,v0,alpha,tols);
Xf = X(end);
Yf = Y(end);
THf = TH(end);
t1f = t1*(J{end}.');
t2f = t2*(J{end}.');
nf = cross(t1f,t2f);
nf = nf/norm(nf);

%% objective function: det(d r/ d(T,\beta))
    function val = objective(g0)
        t0 = T0 + ep0*cos(g0);
        b0 = beta0 + ep0*sin(g0);
        q0 = q_eq + ep*(cos(b0)*t1 + sin(b0)*t2); % initial condition
        [X_,Y_,TH_,J_,~,T_] = swimmerTangentU_Jt(q0(1),q0(2),q0(3),[0 t0],flow,v0,alpha,tols);
        dXdT = flow.Ux([X_(end),Y_(end)],T_(end)) + v0*cos(TH_(end));
        dYdT = flow.Uy([X_(end),Y_(end)],T_(end)) + v0*sin(TH_(end));
        dq0dBeta = ep*(-sin(b0)*t1.' + cos(b0)*t2.');
        dRdBeta = J_{end}(1:2,1:3)*dq0dBeta;
        dRmat = [dXdT,dYdT; dRdBeta.'];
        val = det(dRmat);
%         val = det(dRmat)/norm(dRdBeta);
                
    end
end