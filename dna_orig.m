function dYdt = dna_orig(t,Y,p)  % TODO - Write the function declaration. Name the function SIRmodel
G2T = Y(1);
RT = Y(2);
mass = Y(3);
G2R = (2*RT*G2T)/(RT + G2T + .001 + sqrt((RT + G2T + .001)^2 - 4*RT*G2T));

% parameters in equation
k1 = p(1);
k2 = p(2);
k3 = p(3);
k4 = p(4);
kp = p(5);
Kmp = p(6);
k2_ = p(7); % represents k2'

%Mass
dmassdt = 0.00495 * mass;

%A6
dG2Tdt= k1 - k2*G2T - k2_*G2R; 

%A7
dRTdt= k3 - k4*RT - (kp*(RT-G2R)*(G2T-G2R)*mass)/(Kmp + RT - G2R);

% TODO - Create output column vector dYdt
dYdt = [dG2Tdt; dRTdt; dmassdt];

end