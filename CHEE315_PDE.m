function [F, x, t, u, constants] = CHEE315_PDE
%PDEX1  Example 1 for PDEPE
%   This is a simple example with an analytical solution that illustrates
%   the straightforward formulation, computation, and plotting of the solution
%   of a single PDE.
%
%   In the form expected by PDEPE, the single PDE is
%
%   [pi^2] .*  D_ [u] =  D_ [   Du/Dx   ] +      [0]
%              Dt        Dx
%    ----         ---       -------------    -------------
%     c            u        f(x,t,u,Du/Dx)   s(x,t,u,Du/Dx)
%
%   The equation is to hold on an interval 0 <= x <= 1 for times t >= 0.
%   There is an analytical solution u(x,t) = exp(-t)*sin(pi*x). Boundary
%   conditions are specified so that it will be the solution of the initial-
%   boundary value problem.  Obviously the initial values are sin(pi*x).
%   Two kinds of boundary conditions are chosen so as to show how they
%   appear in the form expected by the solver.  The left bc is u(0,t) = 0:
%
%     [u]    +     [0] .* [       Du/Dx          ] = [0]
%
%     ---          ---    ------------------------   ---
%   p(0,t,u)      q(0,t)        f(0,t,u,Du/Dx)        0
%
%   The right bc is taken to be pi*exp(-t) + Du/Dx(1,t) = 0:
%
%   [ pi*exp(-t) ] + [1] .* [   Du/Dx    ] = [0]
%
%   --------------   ---    --------------   ---
%      p(1,t,u)     q(1,t)   f(1,t,u,Du/Dx)   0
%
%   The problem is coded in subfunctions PDEX1PDE, PDEX1IC, and PDEX1BC.
%
%   See also PDEPE, FUNCTION_HANDLE.

%   Lawrence F. Shampine and Jacek Kierzenka
%   Copyright 1984-2014 The MathWorks, Inc.



%Important variables to define
%To put a number on my variable I decided to imagine the plain wall made of
%stainless steel cubes that have 5 cm sides and each generates 400 Watts 
%due to a current that is flown through them. The width of the wall is L =
%0.05 meters
global qgen rho cp k h Tinf

%Heat generation rate of one cube
V = 0.05^3; % m^3   
q = 400; % Watts 
qgen = q/V; %W/m^3

%AISI 304 Stainless steel properties at 300K (I will keep these values 
%constant altough they change with temperature)
rho = 7900; %kg/m^3
cp  = 477; %J/(kg*K)
k   = 14.9; %W/(m*K)

%Forced convection coefficient and Temperature
Tinf = 40; %degrees celcius
h = 1000; %W/(K*m^2) a low forced convection coefficient for a liquid

constants = struct('qgen',qgen,'rho',rho,'cp',cp,'k',k,'Tinf',Tinf,'h',h);

% The plain wall has width 0.05 meters so L = 0.05 m
% x is a vector that defines the mesh points of my space. So 101 mesh
% points from 0 to 0.05 m
wallL = 0.05;
x = linspace(0,wallL,101);
% t is a vector that defines the time points of the simulation. So 201
% points from 0 toi 500 seconds
tpoints = 201;
t = linspace(0,500,tpoints);

m = 0;
%This solves my PDE
sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);
% Extract the first solution component as u.  This is not necessary
% for a single equation, but makes a point about the form of the output.
u = sol(:,:,1);

%loop over all time points of the simulation. 
loops = tpoints;

%There is need to create a struct variable that will hold the data that 
%makes up the figure for each time point.
F(loops) = struct('cdata',[],'colormap',[]); 

%Get the axis limits of the last time point. 
figure;
plot(x,u(end,:));
xl = xlim;
yl = ylim;

%Run the loop to create all the plots
for j = 1:loops
    tval = j;
    plot(x,u(tval,:),'LineWidth',2);
    hold on
    plot(wallL,Tinf,'*r','MarkerSize',10);
    text(wallL-0.005,Tinf, 'T\infty \rightarrow', 'FontSize',13);
    hold off
    set(gca,'FontSize',12');
    xlim(xl);
    ylim(yl);
    xlabel('Wall length (m)','FontSize',12);
    ylabel('Temperature (^{o}C)','FontSize',12);
    title('Temperature distribution in a 1D wall with heat generation','FontSize',12);
    l = legend(sprintf('Temperature at t = %4.1f, seconds',t(tval)));
    l.FontSize = 12;
    F(j) = getframe(gcf);
end


% --------------------------------------------------------------------------
%This subfunction defines the PDE that is necessary to solve
function [c,f,s] = pdex1pde(x,t,u,DuDx)
global rho cp k qgen;
c = rho*cp/k;
f = DuDx;
s = qgen/k;

% --------------------------------------------------------------------------
%This function defines the initial condition (in our case the plain wall is
%at 0 degrees celcius)
function u0 = pdex1ic(x)
u0 = 0;

% --------------------------------------------------------------------------
%This function defines the boundary conditions
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
global Tinf k h;
pl = ul;
ql = 0;
pr = ur-Tinf;
qr = k/h;


