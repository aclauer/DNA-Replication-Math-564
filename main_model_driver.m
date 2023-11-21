% Script to solve and plot main model
clear all

tspan = linspace(0, 100, 100);

% Initial conditions
y0 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];

% Solve via ode45
[tsol, ysol] = ode45(@main_model, tspan, y0);

plot(tsol, ysol);

xlabel('Time');
ylabel('Concentration');

% TODO: Create legend with selected species