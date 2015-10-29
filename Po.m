function [Po] = Po(X, xa)
deltaE  = 65e-21;       % 65 zeptojoules
d       = 7e-9;         % 7 nm
gamma   = 0.14; 
xc      = 12e-9;        % 12 nm
k       = 1.38065e-23;  % J/K, Boltzmann
Temp    = 293.2;        % Kelvin

Po = 1 / [ 1 + exp((deltaE - Kgs * d * (gamma * X - xa + xc - d/2)) / (k*T))];

