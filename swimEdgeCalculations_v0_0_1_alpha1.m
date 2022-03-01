clear
v0 = 0.1;
alpha = 1;
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle')
load('flows/vortexTimeIndep','flow')
q_eq = [0,asin(v0)/(2*pi),pi/2];
% t1 = [0,0,1];
% t2 = [1,0,0];
t1 = [0.843038096720760,0,-0.537853853270047]; % weak SwIM direction
t2 = [0.016158485134762,0,-0.999869443156530]; % strong SwIM direction


%% set initial guesses based on my numerical investigations
Ti = zeros(1,4);
beta0 = zeros(1,4); 
% first two frontiers
load('swimedge/v0_0_1_alpha_1_Tf_3_beta_range_first_two_frontiers.mat','beta_guess_refined','T_change')
Ti(1:2) = T_change(1:2);
beta0(1:2) = beta_guess_refined(1:2);
% third and fourth frontiers
load('swimedge/v0_0_1_alpha_1_Tf_3_beta_range_third_frontier.mat','beta_guess_refined','T_change')
Ti(3:4) = T_change(2:3);
beta0(3:4) = beta_guess_refined(2:3);


ep = 1e-3;
tols = 1e-9*[1,1];

dT = 0.1; % initial change in integration time
Tf = 4; % final integration time

% output variables
betac = cell(1,length(Ti));
Xfc = cell(1,length(Ti));
Yfc = cell(1,length(Ti));
THfc = cell(1,length(Ti));
nfc = cell(1,length(Ti));
fvalc = cell(1,length(Ti));
exitFlagc = cell(1,length(Ti));
outputc = cell(1,length(Ti));
dobjc = cell(1,length(Ti));
Tc = cell(1,length(Ti));

for i = 1:length(Ti)
    [betac{i},Xfc{i},Yfc{i},THfc{i},nfc{i},fvalc{i},exitFlagc{i},outputc{i},dobjc{i},Tc{i}] = swimEdgeContinuationXY(v0,alpha,flow,q_eq,t1,t2,beta0(i),Ti(i),Tf);
end

% plot results
figure
subplot(1,2,1)
% plot SwIM edge determined by picking innermost points of patch of stable
% manifold
load('barrier_analysis/swim2D_eSwIM','XeSwim','YeSwim')
% transform from stable to unstable manifold by reflection about y = x
Xu = YeSwim;
Yu = -XeSwim; % transform 
font = 22;
plot(Xu,Yu,'k:','LineWidth',4)
% plot folds calculated from continuation
for i = 1:length(Ti)
    hold on
    plotSwimmerTrajectory(Xfc{i},Yfc{i},atan2(nfc{i}(:,2),nfc{i}(:,1)),v0,alpha,flow,1,'o-')
end
set(gca,'FontSize',font)
ylim([-0.05,0.55])
xlim([-0.05,0.55])


% Plot just innermost pieces of folds
subplot(1,2,2)
% plot SwIM edge determined by picking innermost points of patch of stable
% manifold
plot(Xu,Yu,'k:','LineWidth',4)
% plot pieces of innermost pieces of each fold 
for i = 1:length(Ti)
    hold on
    switch i
        case 1
            inds = 1:34;
        case 2
            inds = 10:17;
        case 3
            inds = 15:29;
        case 4
            inds = 10:length(Xfc{4});
    end
    plotSwimmerTrajectory(Xfc{i}(inds),Yfc{i}(inds),atan2(nfc{i}(inds,2),nfc{i}(inds,1)),v0,alpha,flow,1)
end
set(gca,'FontSize',font)
ylim([-0.05,0.55])
xlim([-0.05,0.55])

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

%% calculate trajectories to illustrate blocking behavior
n=1e4;
tf=3;
% initialize trajectories on line y = 0, with theta in (0,pi) so everyone
% swimming up
X0 = 1*rand([n,1]) - 0.5;
Y0 = zeros([n,1]);
TH0 = pi*rand([n,1]);
tols = [1e-9,1e-9];
[X,Y,TH,T] = swimmerTrajectoryU(X0,Y0,TH0,[0,tf],flow,v0,alpha,tols);

%% plot trajectories illustrating blocking behavior, in two panels with different cutoff times
% collate folds into single curve: the SwIM edge
    % plot computed SwIM edge
Xswim = []; Yswim = []; THswim = [];
    for i = 1:length(Ti)
        switch i
            case 1
                inds = 1:34;
            case 2
                inds = 10:17;
            case 3
                inds = 15:29;
            case 4
                inds = 10:length(Xfc{4});
        end
        Xswim = [Xswim; Xfc{i}(inds)];
        Yswim = [Yswim; Yfc{i}(inds)];
        THswim = [THswim;atan2(nfc{i}(inds,2),nfc{i}(inds,1))];
    end

figure
for iter = 1:2
    subplot(1,2,iter)
    % plot BIM
    plotSwimmerTrajectory(Xs,Ys,THs,v0,alpha,flow,1,'r:')
    hold on
    % plot reflection of BIM
    plotSwimmerTrajectory(-Xs,Ys,pi-THs,v0,alpha,flow,0,'r:')
    hold on
    % plot computed SwIM edge
    plotSwimmerTrajectory(Xswim,Yswim,THswim,v0,alpha,flow,0,'r',0)
    hold on
    % plot reflection
    plotSwimmerTrajectory(-Xswim,Yswim,pi-THswim,v0,alpha,flow,0,'r',0)
    hold on
    
    if iter == 1
        tcut = 1.66;
    else
        tcut = tf;
    end
    % plot trajectories
    plot(X(T <= tcut,:),Y(T <= tcut,:),'g')
    if iter == 1
        label = ['(a) $t \leq ',num2str(tcut),'$'];
    else
        label = ['(b) $t \leq ',num2str(tcut),'$'];
    end
    text(-0.9,0.6,label,'Interpreter','latex','FontSize',font)
    
    ylim([-0.1,0.6])
    xlim([-1,1])
    set(gca,'YTick',[0,0.5],'XGrid','on','YGrid','on','GridAlpha',0.6,'FontSize',font)
end

% ylim([-0.05,0.55])
% xlim([-0.05,0.55])
