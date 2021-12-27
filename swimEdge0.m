% SWIMEDGE0 searches for the initial condition of a point hitting the SwIM
% edge at time T. Initial condition of form q0 = q_eq + eps*(cos(beta)*t1 +
% sin(beta)*t2)
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
% fval, exitFlag, output - fsolve outputs
% dobj - derivative of the objective function at the solution
function [beta,Xf,Yf,THf,nf,fval,exitFlag,output,dobj] = swimEdge0(beta0,T,ep,q_eq,t1,t2,v0,alpha,flow,tols)
%% compute intial condition of trajectory going to SwIM edge in time T
opt = optimoptions('fsolve','Display','iter-detailed');
[beta,fval,exitFlag,output,dobj] = fsolve(@objective,beta0,opt);

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
        [X_,Y_,TH_,J_,~,~] = swimmerTangentU_Jt(q0(1),q0(2),q0(3),[0 T],flow,v0,alpha,tols);
        t1T = t1*(J_{end}.');
        t2T = t2*(J_{end}.');
        nT = cross(t1T,t2T);
        val = nT(3);%./norm(nT); % theta component of t1T X t2T, unnormalized at present
    end
end