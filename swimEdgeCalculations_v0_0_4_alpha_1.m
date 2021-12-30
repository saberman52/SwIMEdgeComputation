clear

%% compute SwIM edge
alpha = 1;
v0 = 0.4;
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle')
load('flows/vortexTimeIndep.mat','flow') % fluid velocity field
q_eq = [0,asin(v0)/(2*pi),pi/2]; % SFP
% compute eigenvectors of SFP
[V,D] = eig(flow.A(q_eq(1),q_eq(2),q_eq(3),v0,alpha,0));
[~,indSort] = sort(diag(D),'descend');
t1 = V(:,indSort(2)).'; % weak SwIM direction
t2 = V(:,indSort(1)).'; % strong SwIM
% make sure eigenvectors point to the right (in +x direction)
if t1(1) < 0
    t1 = -t1;
end
if t2(1) < 0
    t2 = -t2;
end

beta0 = 0; % initial guess for fold: fully along weak SwIM direction
Ti = 0.2; % initial integration time
Tf = 2; % max integration time
dT = 0.025; % time step along fold

% compute SwIM edge
[betac,Xfc,Yfc,THfc,nfc,fvalc,exitFlagc,outputc,dobjc,Tc] = swimEdgeContinuation(v0,alpha,flow,q_eq,t1,t2,beta0,Ti,Tf,dT);

%% compute BIM for same fixed point
aBim = -1;
[V,D] = eig(flow.A(q_eq(1),q_eq(2),q_eq(3),v0,aBim,0));
[~,indMax] = max(diag(D));
bimDir = V(:,indMax);
epsilon = 1e-3;
if bimDir(1) < 0
    bimDir = -bimDir;
end
qS0 = q_eq.' + epsilon*bimDir;
tspan = [0 Tf];
tols = 1e-9*[1,1];
[Xs,Ys,THs,Ts] = swimmerTrajectoryU(qS0(1),qS0(2),qS0(3),tspan,flow,v0,aBim,tols);

%% plot SwIM edge, BIM, and boundary of slow zone around vortex center (where |u| = v0) [panel (a)]
figure
% plot SwIM edge with manifold normal vector
plotSwimmerTrajectory(Xfc,Yfc,atan2(nfc(:,2),nfc(:,1)),v0,alpha,flow,1,'r')
hold on
% plot BIM
plotSwimmerTrajectory(Xs,Ys,THs,v0,aBim,flow,0,'r:',0)
set(gca,'FontSize',24)
title('')

%% calculate trajectories to illustrate blocking behavior
n=1e3;
tf=0.7;
% initialize trajectories on line y = 0, with theta in (0,pi) so everyone
% swimming up
X0 = 1*rand([n,1]) - 0.5;
Y0 = zeros([n,1]);
TH0 = pi*rand([n,1]);
tols = [1e-9,1e-9];
[X,Y,TH,T] = swimmerTrajectoryU(X0,Y0,TH0,[0,tf],flow,v0,alpha,tols);

%% plot trajectories illustrating blocking behavior [panel (b)]
figure
% plot computed SwIM edge
plotSwimmerTrajectory(Xfc,Yfc,atan2(nfc(:,2),nfc(:,1)),v0,alpha,flow,1,'r')
hold on
% plot reflection of this SwIM edge
plotSwimmerTrajectory(-Xfc,Yfc,pi - atan2(nfc(:,2),nfc(:,1)),v0,alpha,flow,0,'r')
hold on
% plot trajectories from monte carlo simulation
plot(X,Y,'g')
ylim
set(gca,'FontSize',24)
title('')
ylim([-0.1,0.6])
set(gca,'YTick',[0,0.5],'XGrid','on','YGrid','on','GridAlpha',0.6)

%% calculate a patch of 2D SwIM
N = 100;
ep = 1e-3; % this is the same value of ep as in the SwIM edge continuation function
beta0 = linspace(-1,0.2,N).'; % initial angle parameter
Q0 = repmat(q_eq,N,1) + ep*(cos(beta0)*t1 + sin(beta0)*t2);
[X,Y,TH,T] = swimmerTrajectoryU(Q0(:,1),Q0(:,2),Q0(:,3),[0,2],flow,v0,alpha,1e-9*[1,1]);

%% plot the 2D SwIM in xy\theta phase space
figure
surf(X,Y,TH,'FaceAlpha',0.35,'LineStyle','none')
hold on
plot3(X,Y,TH,'r')
% plot SwIM edge in phase space
hold on
plot3(Xfc,Yfc,THfc,'ko-','LineWidth',2)
set(gca,'FontSize',24)
xlim([0,0.4])
ylim([0,0.4])
zlim([0,6])
colorbar
xlabel('x')
ylabel('y')
zlabel('\theta')
% view([18.52 33.53])