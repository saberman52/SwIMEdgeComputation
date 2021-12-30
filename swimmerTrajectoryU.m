% SWIMMERTRAJECTORYU integrates the trajectories of swimmers in a 2D flow
% field. The swimmers are noninteracting.
% Input:
% [X0,Y0] - vectors of initial positions of swimmers
% [TH0] - vector of initial angles of swimmers
% tspan - [t0 tf], initial and final time of integration
% flow - struct of functions specifying flow field (not necessarily incompressible):
%   flow.Ux(r,t) - x component of flow velocity
%   flow.Uy(r,t) - y component of flow velocity
%   flow.Uxx(r,t) - x derivative of Ux
%   flow.Uxy(r,t) - y derivative of Ux
%   flow.Uyx(r,t) - x derivative of Uy
%   flow.Uyy(r,t) - y derivative of Uy
%   r - matrix of (x,y) pairs, where first column is x and second column is
%   y
% tols - [reltol,abstol] tolerances for integrator (optional)
% Output:
% [X,Y] - matrices of trajectories of each swimmer
% TH0 - matrix of angle trajectory of each swimmer
% T - vector of time steps
function [X,Y,TH,T] = swimmerTrajectoryU(X0,Y0,TH0,tspan,flow,v0,alpha,tols)
    if nargin < 8
        tols = [1e-8,1e-6]; % to be verified
    end
    options = odeset('RelTol',tols(1),'AbsTol',tols(2));
    Z0 = zeros(length(X0)*3,1);
    Z0(1:3:end) = X0;
    Z0(2:3:end) = Y0;
    Z0(3:3:end) = TH0;
    [T,Z] = ode45(@V,tspan,Z0,options);
    X = Z(:,1:3:end);
    Y = Z(:,2:3:end);
    TH = Z(:,3:3:end);
    
    function dZ = V(t,Z)
        R = [Z(1:3:end),Z(2:3:end)];
        dZ = zeros(length(Z),1);
        dZ(1:3:end) = flow.Ux(R,t) + v0*cos(Z(3:3:end)); % dX
        dZ(2:3:end) = flow.Uy(R,t) + v0*sin(Z(3:3:end)); % dY
        dZ(3:3:end) = -0.5*(flow.Uxy(R,t) - flow.Uyx(R,t)) ...
            + alpha*(0.5*(flow.Uxy(R,t)+flow.Uyx(R,t)).*cos(2*Z(3:3:end)) + flow.Uyy(R,t).*sin(Z(3:3:end)).*cos(Z(3:3:end)) ...
            - flow.Uxx(R,t).*sin(Z(3:3:end)).*cos(Z(3:3:end)));% dTH
    end
end