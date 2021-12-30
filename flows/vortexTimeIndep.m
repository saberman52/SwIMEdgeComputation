clear
flow = struct;
flow.Ux = @(R,t) sin(2*pi*R(:,1)).*cos(2*pi*R(:,2));
flow.Uy = @(R,t) -cos(2*pi*R(:,1)).*sin(2*pi*R(:,2));
flow.Uxx = @(R,t) 2*pi*cos(2*pi*R(:,1)).*cos(2*pi*R(:,2));
flow.Uxy = @(R,t) -2*pi*sin(2*pi*R(:,1)).*sin(2*pi*R(:,2));
flow.Uyx = @(R,t) 2*pi*sin(2*pi*R(:,1)).*sin(2*pi*R(:,2));
flow.Uyy = @(R,t) -2*pi*cos(2*pi*R(:,1)).*cos(2*pi*R(:,2));

flow.A = @(x,y,th,v0,alpha,t) [2*pi*cos(2*pi*x)*cos(2*pi*y), -2*pi*sin(2*pi*x)*sin(2*pi*y), -v0*sin(th);...
                                2*pi*sin(2*pi*x)*sin(2*pi*y), -2*pi*cos(2*pi*x)*cos(2*pi*y), v0*cos(th);...
4*pi^2*(cos(2*pi*x)*sin(2*pi*y) + alpha*sin(2*pi*x)*cos(2*pi*y)*sin(2*th)), 4*pi^2*(sin(2*pi*x)*cos(2*pi*y) + alpha*cos(2*pi*x)*sin(2*pi*y)*sin(2*th)),...
-4*pi*alpha*cos(2*pi*x)*cos(2*pi*y)*cos(2*th)]; % only defined for individual initial conditions, not vector

% velocity gradient matrix for swimmer equations in 4D phase space
flow.A4 = @(x,y,nx,ny,v0,alpha,t) [2*pi*cos(2*pi*x)*cos(2*pi*y), -2*pi*sin(2*pi*x)*sin(2*pi*y), v0/sqrt(nx^2+ny^2)*(1 - nx^2/(nx^2+ny^2)), -v0*nx*ny/(nx^2+ny^2)^1.5;...
                                    2*pi*sin(2*pi*x)*sin(2*pi*y), -2*pi*cos(2*pi*x)*cos(2*pi*y),  -v0*nx*ny/(nx^2+ny^2)^1.5, v0/sqrt(nx^2+ny^2)*(1 - ny^2/(nx^2+ny^2));...
4*pi^2*(-alpha*sin(2*pi*x)*cos(2*pi*y)*nx-cos(2*pi*x)*sin(2*pi*y)*ny), 4*pi^2*(-alpha*cos(2*pi*x)*sin(2*pi*y)*nx - sin(2*pi*x)*cos(2*pi*y)*ny), 2*pi*alpha*cos(2*pi*x)*cos(2*pi*y), -2*pi*sin(2*pi*x)*sin(2*pi*y);...
4*pi^2*(cos(2*pi*x)*sin(2*pi*y)*nx+alpha*sin(2*pi*x)*cos(2*pi*y)*ny), 4*pi^2*(sin(2*pi*x)*cos(2*pi*y)*nx + alpha*cos(2*pi*x)*sin(2*pi*y)*ny), 2*pi*sin(2*pi*x)*sin(2*pi*y), -2*pi*alpha*cos(2*pi*x)*cos(2*pi*y)];

save('vortexTimeIndep')