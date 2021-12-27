% SWIMEDGE0_det_fixedT searches for the initial condition of a point hitting the SwIM
% edge at time T. Initial condition of form q0 = q_eq + eps*(cos(beta)*t1 +
% sin(beta)*t2).
% Objective function: determinant of dr/d(T,beta)
% Inputs:
% beta0 - initial guess for initial condition parameter
% T - integration time (fixed)
% ep - initial step size away from equilibrium point
% q_eq - equilibrium point (row vector)
% t1 - one of the tangent vectors (row vector)
% t2 - other tangent vector
% v0 - swim speed
% alpha - swimmer shape parameter
% flow - fluid velocity field
% tols - numerical integration tolerances
%
% Outputs:
% beta - final initial condition parameter
% Xf - final x coordinate on SwIM edge
% Yf - final y coordinate
% THf - final theta coordinate
% nf - final normal vector to manifold (normalized)
% fval,exitFlag,output - outputs of fsolve/fzero
% dobj - derivative of the objective function at the solution
function [beta,Xf,Yf,THf,nf,fval,exitFlag,output,dobj] = swimEdge0_det_fixedT(beta0,T,ep,q_eq,t1,t2,v0,alpha,flow,tols)
%% compute intial condition of trajectory going to SwIM edge in time T
%% fsolve
% opt = optimoptions('fsolve','Display','iter-detailed','OptimalityTolerance',1e-8);
% [beta,fval,exitFlag,output,dobj] = fsolve(@objective,beta0,opt);
%% fzero
opt = optimset('Display','iter');
[beta,fval,exitFlag,output] = fzero(@objective,beta0,opt);
dobj = 0;

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

%% objective function
    function val = objective(b0)
        q0 = q_eq + ep*(cos(b0)*t1 + sin(b0)*t2); % initial condition
        [X_,Y_,TH_,J_,~,T_] = swimmerTangentU_Jt(q0(1),q0(2),q0(3),[0 T],flow,v0,alpha,tols);
        dXdT = flow.Ux([X_(end),Y_(end)],T_(end)) + v0*cos(TH_(end));
        dYdT = flow.Uy([X_(end),Y_(end)],T_(end)) + v0*sin(TH_(end));
        dq0dBeta = ep*(-sin(b0)*t1.' + cos(b0)*t2.');
        dRdBeta = J_{end}(1:2,1:3)*dq0dBeta;
        dRmat = [dXdT,dYdT; dRdBeta.'];
        val = det(dRmat);
%         val = det(dRmat)/norm(dRdBeta);
    end
end