clear
% load flow
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle')
load('flows/vortexTimeIndep.mat','flow')

v0 = 0.4;
q_eq = [0,asin(v0)/(2*pi),pi/2]; % SFP
% SwIM tangent plane at SFP
t1 = [0,0,1];
t2 = [1,0,0];

% parameters controlling SwIM edge continuation
arcmax = 4;
ep0 = 0.05;

%% alpha = 0.1
alpha = 0.1;
% primary fold
T0 = 0.2;
beta0 = 3;
dT = 0.05;
[gammac,betac,Xfc,Yfc,THfc,nfc,fvalc,exitFlagc,outputc,dobjc,Tc] = swimEdgeContinuationCircDet(v0,alpha,flow,q_eq,t1,t2,beta0,T0,arcmax,dT,ep0);
save('swimedge/lowerfold_v0_4_alpha0_1.mat')

% secondary fold
T0 = 1;
beta0 = 3.05;
dT = -0.05;
[gammac,betac,Xfc,Yfc,THfc,nfc,fvalc,exitFlagc,outputc,dobjc,Tc] = swimEdgeContinuationCircDet(v0,alpha,flow,q_eq,t1,t2,beta0,T0,arcmax,dT,ep0);
save('swimedge/upperfold_v0_4_alpha0_1.mat')

%% alpha = 0.3
alpha = 0.3;
% primary fold
T0 = 0.2;
beta0 = 3;
dT = 0.05;
[gammac,betac,Xfc,Yfc,THfc,nfc,fvalc,exitFlagc,outputc,dobjc,Tc] = swimEdgeContinuationCircDet(v0,alpha,flow,q_eq,t1,t2,beta0,T0,arcmax,dT,ep0);
save('swimedge/lowerfold_v0_4_alpha0_3.mat')

% secondary fold
T0 = 1;
beta0 = 2;
dT = -0.05;
[gammac,betac,Xfc,Yfc,THfc,nfc,fvalc,exitFlagc,outputc,dobjc,Tc] = swimEdgeContinuationCircDet(v0,alpha,flow,q_eq,t1,t2,beta0,T0,arcmax,dT,ep0);
save('swimedge/upperfold_v0_4_alpha0_3.mat')

%% plot
plot_multiple_swim_edges
