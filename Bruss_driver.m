clear all
%close all

% First solve the ode system

%%%%%%%%%%%%%%%%%%%%%%
% INITIAL CONDITION
Y0 = [0.5;1.5];
%%%%%%%%%%%%%%%%%%%%%%

%timespan
tRange = [0 100];

% parameters in Brusselator
k1=1.3;
k2=1;
k3=1/2;
k4=1;

p=[k1,k2,k3,k4];

% call the ODE solver ode15s instead of ode45
% to send parameters to the ode solver, use the following command:
[tSol,YSol] = ode15s(@(tSol,YSol)bruss(tSol,YSol,p),tRange,Y0);


% plot solutions in time
figure(1)
clf
plot(tSol,YSol,'LineWidth',2)
xlabel('Time')
ylabel('Concentration')
legend('S_1','S_2')
set(gca,'FontSize',18)
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the phase portrait
%%%%%%%%%%%%%%%%%%%%%%%%%

% generate grid for phase plane [0 M] X [0 N]
M=5;
N=5;
[xx,yy]=meshgrid(0:0.25:M,0:0.25:N);


% define the right hand sides
ds1dt=k1-k2*xx+k3*xx.^2.*yy - k4*xx;
ds2dt=k2*xx-k3*xx.^2.*yy;

% L = scale factor to make all of the arrows the same length
L = sqrt(ds1dt.^2 + ds2dt.^2);

% plot the direction field
figure(2)
clf % clear the current figure
quiver(xx,yy,ds1dt./L,ds2dt./L,0.5,'LineWidth',1.5); % The final argument is a scaling factor
hold on

% define the nullclines
ncfun1 = @(x,y) k1-k2*x+k3*x.^2.*y - k4*x;
ncfun2 = @(x,y) k2*x-k3*x.^2.*y;
% plot the nullclines
fimplicit(@(x,y) ncfun1(x,y), [0 M 0 N],'LineWidth',3)
fimplicit(@(x,y) ncfun2(x,y), [0 M 0 N],'LineWidth',3);

% plot the trajectories in blue
% these are just plotting s1 vs s2 using solutions from above
plot(YSol(:,1),YSol(:,2),'b','LineWidth',2)
axis tight