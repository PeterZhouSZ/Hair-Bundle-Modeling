function [ CaOUT,terfc ] = CaOUTf(tspan,onset,offset,t)
I = 30e-9;
r = 2e-6;
bathCa = 1e-6;
Ca = zeros(length(tspan));
terfc = (tspan<=onset).*(onset+1) + ((tspan>onset) & (tspan<offset)).*tspan + (tspan>=offset).*(onset+1);

%Ca = (tspan<=onset).* bathCa + ((tspan>onset) & (tspan<offset)).*(bathCa + (((0.12 * I) / (4 * pi * 8e-10 * 1.602e-19 * r))) * erfc(r./(2*(sqrt(8e-10*(terfc-onset))))) / 6.022e23) + (tspan>=offset).* (bathCa);
CaOUT = (tspan<=onset).* bathCa + ((tspan>onset) & (tspan<offset)).*(bathCa + (((0.12 * I) / (4 * pi * 8e-10 * 1.602e-19 * r))) * erfc(r./(2*(sqrt(8e-10*(terfc-onset))))) / 6.022e23) + (tspan>=offset).* bathCa;
%CaOUTinterp = interp1(tspan,CaOUT,t)

%tinterp = [0:0.0001:0.2];
%[CaOUTInterp] = interp1(tspan,CaOUT,tinterp);
%plot(tspan,Ca,'o',tinterp,CaOUTInterp)

%[Ca] = 0;
%if t<=onset; Ca = bathCa;
%elseif (t>onset) & (t<offset);
%Ca = bathCa + (((0.12 * I) / (4 * pi * 8e-10 * 1.602e-19 * r))) * erfc(r./(2*(sqrt(8e-10*(t-onset))))) / 6.022e23;
%elseif t>=offset;
%Ca(region3) = bathCa;
%end
