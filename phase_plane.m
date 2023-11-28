function dYdt = phase_plane(t,Y,p)  % TODO - Write the function declaration. Name the function SIRmodel
% TODO - Extract S1, S2
G2T = Y(1);
RT = Y(2);


% TODO - Define the constants
% parameters in Brusselator
k1=p(1);
k2=p(2);
k2_=p(3);
k3=p(4);
k4=p(5);
kp=p(6);
Kmp=p(7);

% G2R = -1;
% G2R as defined by eqn. A5
G2R = (2 * RT * G2T) / (RT + G2T + (((RT - G2R) * (G2T - G2R)) / G2R) + sqrt((RT + G2T + (((RT - G2R) * (G2T - G2R)) / G2R))^2 - 4 * RT * G2T));

% lambda as defined by eqn. A4
l = ((RT - G2R) * (G2T - G2R)) / G2R;

% TODO: Define mass correctly
mass = 1;

dG2Tdt=k1 - k2 * G2T - k2_ * G2R;
dRTdt=k3 - k4 * RT - (kp * (RT - G2R) * (G2T - G2R) * mass)/(Kmp + RT - G2R);


% TODO - Create output column vector dYdt
dYdt = [dG2Tdt; dRTdt];
end