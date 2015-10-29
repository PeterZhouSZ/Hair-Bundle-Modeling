clc
I = 30e-9;  % in Amps
r = 1e-6;   % in meters

Cinf = ((0.12 * I) / (4 * pi * 8e-10 * 1.602e-19 * r));
Cinf = Cinf / 6.022e23 % answer is in millimolar; if desire M, divide by 1000

t = ((r/erfcinv((4 * pi * 8e-10 * 1.602e-19 * r * 0.95 * Cinf) / (0.12 * I)))^2) / (4 * 8e-10) * 1000 % in milliseconds, remove the 1000 if want sec


%%
clf

for r = [1e-6 2e-6 3e-6 4e-4 5e-6]
I = 30e-9;
T = 10;
%r = 5e-6;
% Time vector
t = [0:T/1000:T];
tt = [T:T/1000:T+T/2];

C =    (((0.12 * I) / (4 * pi * 8e-10 * 1.602e-19 * r)) *   erfc(r./(2*(sqrt(8e-10*t)))) / 6.022e23);
Coff = (((0.12 * I) / (4 * pi * 8e-10 * 1.602e-19 * r)) * ( erfc(r./(2*(sqrt(8e-10*T)))) - erfc(r./(2*sqrt((8e-10)*(tt-T)))) )) /6.022e23;

plot(t,C)
hold on
plot(tt,Coff)
hold on
end
