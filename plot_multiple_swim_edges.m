%% graphs for alpha = 0.1
load('swimedge/lowerfold_v0_4_alpha0_1.mat','Xfc','Yfc','q_eq')
Xlower = [q_eq(1); Xfc];
Ylower = [q_eq(2); Yfc];
load('swimedge/upperfold_v0_4_alpha0_1.mat','Xfc','Yfc')
Xupper = Xfc;
Yupper = Yfc;

line = 'LineWidth';
lw = 2;
mkr = 'MarkerSize';
ms = 18;
fs = 22;
int = 'Interpreter';
la = 'latex';
font = 'FontSize';

figure
subplot(1,2,1)
plot(Xlower,Ylower,'r',line,lw)
hold on
plot(Xupper,Yupper,'b',line,lw)
hold on
plot(q_eq(1),q_eq(2),'r.',mkr,ms)
xlim([  -0.0001   0.017382956542260])
ylim([0.065466456963486   0.065615852724427])
xlabel('$x$',int,la)
ylabel('$y$',int,la)
set(gca,'FontSize',fs)
xt = 1e-3;
yt = 0.0656;
text(xt,yt,'(a) $\alpha = 0.1$',font,fs,int,la)


%% graphs for alpha = 0.3
load('swimedge/lowerfold_v0_4_alpha0_3.mat','Xfc','Yfc')
Xlower = [q_eq(1); Xfc];
Ylower = [q_eq(2); Yfc];
load('swimedge/upperfold_v0_4_alpha0_3.mat','Xfc','Yfc')
Xupper = Xfc;
Yupper = Yfc;

subplot(1,2,2)
plot(Xlower,Ylower,'r',line,lw)
hold on
plot(Xupper,Yupper,'b',line,lw)
hold on
plot(q_eq(1),q_eq(2),'r.',mkr,ms)
xlim([ -0.0001    0.017382956542260])
ylim([0.065466456963486   0.065615852724427])
xlabel('$x$',int,la)
ylabel('$y$',int,la)
set(gca,'FontSize',fs)
text(xt,yt,'(b) $\alpha = 0.3$',int,la,font,fs)
