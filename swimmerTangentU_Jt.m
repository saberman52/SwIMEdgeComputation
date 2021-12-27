% SWIMMERTANGENTU_Jt integrates the tangent flow of swimmers in a 2D flow
% field, for a single trajectory
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
%   flow.A(x,y,th,v0,alpha,t) - matrix of partial derivatives of velocity
%   field
%   r - matrix of (x,y) pairs, where first column is x and second column is
%   y
% tols - [reltol,abstol] tolerances for integrator (optional)
% Output:
% [X,Y] - matrices of trajectories of each swimmer
% TH0 - matrix of angle trajectory of each swimmer
% Jt - tangent flow matrix time series as a cell, where each element is a
% matrix
% detJ - time series of determinant of Jacobian along orbit
% T - vector of time steps
function [X,Y,TH,J,detJ,T] = swimmerTangentU_Jt(X0,Y0,TH0,tspan,flow,v0,alpha,tols)
    if nargin < 8
        tols = [1e-8,1e-6]; % to be verified
    end
    options = odeset('RelTol',tols(1),'AbsTol',tols(2));
    Z0 = zeros(12,1);
    Z0(1) = X0;
    Z0(2) = Y0;
    Z0(3) = TH0;
    J0 = eye(3);
    Z0(4:end) = J0(:);
    [T,Z] = ode45(@V,tspan,Z0,options);
    X = Z(:,1);
    Y = Z(:,2);
    TH = Z(:,3);
    Jcell = mat2cell(Z(:,4:end),ones(1,size(Z,1)));
    J = cellfun(@(x) reshape(x,[3,3]),Jcell,'UniformOutput',0);
    detJ = cellfun(@det,J);
    
    
    function dZ = V(t,Z)
        R = [Z(1),Z(2)];
        dZ = zeros(length(Z),1);
        dZ(1) = flow.Ux(R,t) + v0*cos(Z(3)); % dX
        dZ(2) = flow.Uy(R,t) + v0*sin(Z(3)); % dY
        dZ(3) = -0.5*(flow.Uxy(R,t) - flow.Uyx(R,t)) ...
            + alpha*(0.5*(flow.Uxy(R,t)+flow.Uyx(R,t)).*cos(2*Z(3)) + flow.Uyy(R,t).*sin(Z(3)).*cos(Z(3)) ...
            - flow.Uxx(R,t).*sin(Z(3)).*cos(Z(3)));% dTH
        Jt = reshape(Z(4:end),[3,3]);
        dJ = flow.A(Z(1),Z(2),Z(3),v0,alpha,t)*Jt;
        dZ(4:end) = dJ(:);
    end
end