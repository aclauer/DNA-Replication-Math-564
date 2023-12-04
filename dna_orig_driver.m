clear all
%close all

% First solve the ode system

%%%%%%%%%%%%%%%%%%%%%%
% INITIAL CONDITION
Y0 = [0.48;0.25;1]; %not sure what to make these
%%%%%%%%%%%%%%%%%%%%%%

%timespan
tRange = [0 1000];

% parameters in equation
k1=0.015;
k2=0.0075;
k3=0.09375;
k4=0.1875;
kp=3.25;
Kmp=0.001;
k2_=0.05; % represents k2'
p = [k1, k2, k3, k4, kp, Kmp, k2_];
mass=1;

% call the ODE solver ode15s instead of ode45
% to send parameters to the ode solver, use the following command:
%[tSol,YSol] = ode15s(@(tSol,YSol)dna_orig(tSol,YSol,p),tRange,Y0);


% plot solutions in time
%figure(1)
%clf
%plot(tSol,YSol,'LineWidth',2)
%xlabel('Time')
%ylabel('Concentration')
%legend('S_1','S_2')
%set(gca,'FontSize',18)
%grid on

%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the phase portrait
%%%%%%%%%%%%%%%%%%%%%%%%%

% generate grid for phase plane [0 M] X [0 N]
M=0.5;
N=2;
[xx,yy]=meshgrid(linspace(0,M,25),linspace(0,N,25));
% xx represents RT and yy represents G2T

% define the right hand sides
%right hand sides also depend on these functions:
G2R = (2.*xx.*yy)./(xx + yy + .001 + sqrt((xx + yy + .001).^2 - 4.*xx.*yy));

%A6
dG2Tdt= k1 - k2.*yy - k2_.*G2R; 
%A7
dRTdt= k3-k4.*xx-(kp.*(xx-G2R).*(yy-G2R).*mass)./(Kmp+xx-G2R);%./(Kmp+yy-G2R);
%Mass
dmassdt= 0.00495*mass;

% L = scale factor to make all of the arrows the same length
L = sqrt(dG2Tdt.^2 + dRTdt.^2);

% plot the direction field
figure(2)
clf % clear the current figure
quiver(xx, yy, dRTdt./L, dG2Tdt./L, 0.2); % The final argument is a scaling factor
hold on

% define the nullclines
% x represents RT and y represents G2T
g2rfun = @(x,y) (2.*x.*y)./(x + y + .001 + sqrt((x + y + .001).^2 - 4.*x.*y));
ncfun1 = @(x,y) k1 - k2*y - k2_*g2rfun(x,y);
ncfun2 = @(x,y) k3 - k4*x - (kp*(x-g2rfun(x,y)).*(y-g2rfun(x,y))*mass)./(Kmp + x - g2rfun(x,y));
% plot the nullclines
fimplicit(@(x,y) ncfun1(x,y), [0 M 0 N],'LineWidth',3)
fimplicit(@(x,y) ncfun2(x,y), [0 M 0 N],'LineWidth',3);

% plot the trajectories in blue
% these are just plotting s1 vs s2 using solutions from above
%plot(YSol(:,1),YSol(:,1),'b','LineWidth',2)
%axis tight