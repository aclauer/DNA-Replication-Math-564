function dYdt = bruss(t,Y,p)  % TODO - Write the function declaration. Name the function SIRmodel
% TODO - Extract S1, S2
xx = Y(1);
yy = Y(2);


% TODO - Define the constants
% parameters in Brusselator
k1=p(1);
k2=p(2);
k3=p(3);
k4=p(4);

ds1dt=k1-k2*xx+k3*xx.^2.*yy - k4*xx;
ds2dt=k2*xx-k3*xx.^2.*yy;


% TODO - Create output column vector dYdt
dYdt = [ds1dt; ds2dt];
end