clear all
%close all

% First solve the ode system

%%%%%%%%%%%%%%%%%%%%%%
% INITIAL CONDITION
Y0 = [1;1];
%%%%%%%%%%%%%%%%%%%%%%

%timespan
tRange = [0 100];

% parameters in Brusselator
k1=1;
k2=1;
k2_=1;
k3=1;
k4=1;
kp=1;
Kmp=1;

p=[k1,k2,k2_,k3,k4,kp,Kmp];

% call the ODE solver ode15s instead of ode45
% to send parameters to the ode solver, use the following command:
[tSol,YSol] = ode15s(@(tSol,YSol)phase_plane(tSol,YSol,p),tRange,Y0);


% plot solutions in time
figure(1)
clf
plot(tSol,YSol,'LineWidth',2)
xlabel('Time')
ylabel('Concentration')
legend('G2T','RT')
set(gca,'FontSize',18)
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the phase portrait
%%%%%%%%%%%%%%%%%%%%%%%%%

% generate grid for phase plane [0 M] X [0 N]
M=5;
N=5;
[G2T,RT]=meshgrid(0:0.25:M,0:0.25:N);

% Copied from bruss.m
% lambda as defined by eqn. A4
%l = (((RT - G2R) * (G2T - G2R)) / G2R);

% TODO: Define mass correctly
mass = 1;

% G2R = -1;
% G2R as defined by eqn. A5
G2R = (2 * RT * G2T) / (RT + G2T + (((RT - G2R) * (G2T - G2R)) / G2R) + sqrt((RT + G2T + (((RT - G2R) * (G2T - G2R)) / G2R))^2 - 4 * RT * G2T));

dG2Tdt=k1 - k2 * G2T - k2_ * G2R;
dRTdt=k3 - k4 * RT - (kp * (RT - G2R) * (G2T - G2R) * mass)/(Kmp + RT - G2R);

% L = scale factor to make all of the arrows the same length
L = sqrt(dG2Tdt.^2 + dRTdt.^2);

% plot the direction field
figure(2)
clf % clear the current figure
quiver(G2T,RT,dG2Tdt./L,dRTdt./L,0.5,'LineWidth',1.5); % The final argument is a scaling factor
hold on

% define the nullclines
%ncfun1 = @(x,y) k1-k2*x+k3*x.^2.*y - k4*x;
%ncfun2 = @(x,y) k2*x-k3*x.^2.*y;
% plot the nullclines
fimplicit(@(x,y) ncfun1(x,y), [0 M 0 N],'LineWidth',3)
fimplicit(@(x,y) ncfun2(x,y), [0 M 0 N],'LineWidth',3);

% plot the trajectories in blue
% these are just plotting s1 vs s2 using solutions from above
plot(YSol(:,1),YSol(:,2),'b','LineWidth',2)
axis tight