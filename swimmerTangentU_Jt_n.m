% SWIMMERTANGENTU_Jt_n integrates the tangent flow of swimmers in a 2D flow
% field, for n trajectories
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
%   flow.A(x,y,th,v0,alpha,t) - stability matrix of swimmer equations
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
function [X,Y,TH,J,detJ,T] = swimmerTangentU_Jt_n(X0,Y0,TH0,tspan,flow,v0,alpha,tols)
    if nargin < 8
        tols = [1e-8,1e-6]; % to be verified
    end
    options = odeset('RelTol',tols(1),'AbsTol',tols(2));
    n = length(X0);
    nel = 12; % number of dynamical variables
    Z0 = zeros(nel*n,1);
    % initialize trajectory variables
    Z0(1:nel:end) = X0;
    Z0(2:nel:end) = Y0;
    Z0(3:nel:end) = TH0;
    % initialize tangent flow variables
    J0 = eye(3);
    mod_inds = mod(0:(nel*n-1),nel);
    tflow_inds = find(mod_inds >= 3); % tangent flow dynamical variable indices are 3 - 11 mod 12.
    Z0(tflow_inds) = repmat(J0(:),n,1);
    [T,Z] = ode45(@V,tspan,Z0,options);
    X = Z(:,1:nel:end);
    Y = Z(:,2:nel:end);
    TH = Z(:,3:nel:end);
    Jcell = mat2cell(Z(:,tflow_inds),ones(1,size(Z,1)),9*ones(1,n));
    J = cellfun(@(x) reshape(x,[3,3]),Jcell,'UniformOutput',0);
    detJ = cellfun(@det,J);
    
    
    function dZ = V(t,Z)
        R = [Z(1:nel:end),Z(2:nel:end)];
        TH_ = Z(3:nel:end);
        dZ = zeros(length(Z),1);
        % trajectory differential
        dZ(1:nel:end) = flow.Ux(R,t) + v0*cos(TH_); % dX
        dZ(2:nel:end) = flow.Uy(R,t) + v0*sin(TH_); % dY
        dZ(3:nel:end) = -0.5*(flow.Uxy(R,t) - flow.Uyx(R,t)) ...
            + alpha*(0.5*(flow.Uxy(R,t)+flow.Uyx(R,t)).*cos(2*TH_) + flow.Uyy(R,t).*sin(TH_).*cos(TH_) ...
            - flow.Uxx(R,t).*sin(TH_).*cos(TH_));% dTH
        %matrix differential
        Jt = mat2cell(Z(tflow_inds),9*ones(1,n),1);
        Jt = cellfun(@(x) reshape(x,[3,3]),Jt,'UniformOutput',0);
        Afun = @(x,y,z) flow.A(x,y,z,v0,alpha,t);
        At = cellfun(Afun,num2cell(R(:,1)),num2cell(R(:,2)),num2cell(TH_),'UniformOutput',0);
        dJ = cellfun(@(x,y) x*y,At,Jt,'UniformOutput',0);
        dZ(tflow_inds) = cell2mat(cellfun(@(x) x(:),dJ,'UniformOutput',0));
    end
end