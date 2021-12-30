%PLOTSWIMMERTRAJECTORY plots a 2D swimmer trajectory obtained by numerical
%integration along with the swimming fixed points in the first vortex cell.
%Plots in current figure
function plotSwimmerTrajectory(X,Y,TH,v0,alpha,flow,plotEverything,color,quiv)
fig = gcf;
lw = 3;
sz = 0.05;
msize = 32;
grey = 0.6*ones(1,3);
% plot vortex boundaries
if nargin < 7
    plotEverything=1;
end
if plotEverything
    plot([0,0.5,0.5,0,0],[0,0,0.5,0.5,0],'--','LineWidth',2,'Color',grey)
    hold on
end

if plotEverything
    %% plot swimming fixed points
    rho = asin(v0)/(2*pi);
    % unstable front fixed points
    xUn = [0, rho, 0.5, 0.5-rho];%, 0.5]; fixed point in adjacent vortex cell
    yUn = [rho, 0.5, 0.5-rho, 0];%, 0.5+rho];
    thUn = [pi/2,0, -pi/2, pi];%, pi/2];
    % stable front fixed points
    
    xSt = [rho,0,0.5-rho,0.5];%, 0.5]; 
    ySt = [0,0.5-rho,0.5,rho];%, 1-rho];
    thSt = [pi,pi/2,0,-pi/2];%, pi/2];
    plot(xUn,yUn,'r.','MarkerSize',msize)
    quiver(xUn,yUn,sz*cos(thUn),sz*sin(thUn),0,'r','LineWidth',lw)
    plot(xSt,ySt,'b.','MarkerSize',msize)
    quiver(xSt,ySt,sz*cos(thSt),sz*sin(thSt),0,'b','LineWidth',lw)

    %%plot slow zones
    n = 200;
    [xx,yy] = meshgrid(linspace(0,0.5,n),linspace(0,0.5,n));
    flowSpeed = reshape(sqrt(flow.Ux([xx(:),yy(:)],0).^2 + flow.Uy([xx(:),yy(:)],0).^2),size(xx));
    fig2 = figure;
    M = contour(xx,yy,flowSpeed,[v0 v0]);

    close(fig2)
    figure(fig)
    keep = find(M(1,2:end) < 0.5 & M(1,2:end) > 0 & M(2,2:end) < 0.5 & M(2,2:end) > 0);
    plot(M(1,1+keep),M(2,1+keep),'.','Color',grey)

    %% plot fixed point contours for periodic vortex array
    % if alpha < 0
    %     fixedPointContour = cos(2*pi*xx).^2.*cos(2*pi*yy).^2;
    %     fig2 = figure;
    %     M = contour(xx,yy,fixedPointContour,-v0^2/(2*alpha)*[1,1]);
    %     close(fig2)
    %     figure(fig)
    %     keep = find(M(1,2:end) < 0.5 & M(1,2:end) > 0 & M(2,2:end) < 0.5 & M(2,2:end) > 0);
    %     plot(M(1,1+keep),M(2,1+keep),'x','Color',(128/255)*ones(1,3))
    % end
    xlabel('$x$','Interpreter','latex')
    ylabel('$y$','Interpreter','latex')
    set(gca,'FontSize',18)
end
% plot trajectory
if nargin >= 8
    if ischar(color)
        p = plot(X,Y,color,'LineWidth',lw);
    else
        p = plot(X,Y,'Color',color,'LineWidth',lw);
    end
else
    p = plot(X,Y,'LineWidth',lw);
end
color = get(p,'Color');
% determine spacing based on traversed arclength
s = cumsum(sqrt(diff(X).^2 + diff(Y).^2));
% spacing = find(s > 0.025,1);
% quiver(X(1:spacing:end),Y(1:spacing:end),sz*cos(TH(1:spacing:end)),sz*sin(TH(1:spacing:end)),0,'LineWidth',lw,'Color',color)
if nargin < 9
    quiv = 1;
end
if quiv
    spacing = find(diff(floor(s/0.1)))+1;
    % spacing = 1:length(X);
    quiver(X(spacing),Y(spacing),sz*cos(TH(spacing)),sz*sin(TH(spacing)),0,'LineWidth',lw,'Color',color)
end
hold off
if plotEverything
    title(['$v_0 = ' num2str(v0) ',\,\,\alpha = ' num2str(alpha) '$'],'Interpreter','latex')
end
% xlim([-0.02 0.52])
axis equal
% ylim([-0.02,0.52])
end