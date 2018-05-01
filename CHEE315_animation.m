%%
%Run the PDE solver to get the data to create the animations
[F,x,t,u,cst] = CHEE315_PDE;


%%
fig1 = figure(1);

%Fast movie of heat transfer in 1D wall
nfast = [1,1:201];
movie(fig1,F,nfast,12);


%%
fig1 = figure(1);

%Slow movie of the first 60 seconds
nslow = [1,1:25];
movie(fig1,F,nslow,1);

%------------------------------%
%%
%Plot of select Temperature lines
plot(x,u(end,:));
xl = xlim;
yl = ylim;

plot(x,u(1,:),'LineWidth',2);
hold on
plot(x,u(4,:),'LineWidth',2);
plot(x,u(16,:),'LineWidth',2);
plot(x,u(25,:),'LineWidth',2);


set(gca,'FontSize',10');
xlabel('Wall length (m)','FontSize',10);
ylabel('Temperature (^{o}C)','FontSize',10);
title('Temperature distribution in a 1D wall with heat generation','FontSize',12);
t0 = sprintf('Temp at t = %4.1f s',t(1));
t4 = sprintf('Temp at t = %4.1f s',t(4));
t16 = sprintf('Temp at t = %4.1f s',t(16));
t25 = sprintf('Temp at t = %4.1f s',t(25));
t200 = sprintf('Temp at t = %4.1f s',t(200));
l = legend(t0,t4,t16,t25,'Location','northwest');
l.FontSize = 10;


hold off

%%
%----------------------------%
%Plot of heat fluxes 


dTdx_0 = (u(:,2)-u(:,1))/(x(2)-x(1));
qxL = (cst.h).*(u(:,101)-(cst.Tinf));
qx0 = -(cst.k)*dTdx_0;
n=[4 16 25];
fig2 = figure(2);

plot(t,qxL, 'LineWidth',2)
hold on
plot(t,qx0,'LineWidth',2);
plot(t(n),qxL(n),'*k','MarkerSize',8);
plot(t(n),qx0(n),'*k','MarkerSize',8);
grid

hold off
set(gca,'FontSize',10);
xlabel('time (s)','FontSize',10);
ylabel('Heat flux (W/m^{2})','FontSize',10);
th = title('Heat flux through boundaries of 1D wall with heat generation','FontSize',12);

%Fix the title position because it overlaped with the 10^4 text
titlePos = get(th, 'Position');
set(th, 'Position', [titlePos(1), titlePos(2)+6000, titlePos(3)]);

l = legend(sprintf('q"_{x}(L,t)'),sprintf('q"_{x}(0,t)'), 'Location','east');
l.FontSize = 10;


