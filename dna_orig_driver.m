clear all

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

%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the phase portrait
%%%%%%%%%%%%%%%%%%%%%%%%%

% generate grid for phase plane [0 M] X [0 N]
M=0.5;
N=2;
[xx,yy]=meshgrid(linspace(0,M,10),linspace(0,N,15));
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
quiver(xx, yy, dRTdt./L, dG2Tdt./L, 0.1); % The final argument is a scaling factor
hold on

% define the nullclines
% x represents RT and y represents G2T
g2rfun = @(x,y) (2.*x.*y)./(x + y + .001 + sqrt((x + y + .001).^2 - 4.*x.*y));
ncfun1 = @(x,y) k1 - k2*y - k2_*g2rfun(x,y);
ncfun2 = @(x,y) k3 - k4*x - (kp*(x-g2rfun(x,y)).*(y-g2rfun(x,y))*mass)./(Kmp + x - g2rfun(x,y));
% plot the nullclines
fimplicit(@(x,y) ncfun1(x,y), [0 M 0 N],'LineWidth',3)
fimplicit(@(x,y) ncfun2(x,y), [0 M 0 N],'LineWidth',3);