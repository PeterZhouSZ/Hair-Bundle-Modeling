function [ Ca ] = CaINmod4(t,onset,offset)
I = 30e-9;
r = 2e-6;
bathCa = 1e-6;
Ca = zeros(length(t));
terfc = (t<=onset).*(onset+1) + ((t>onset) & (t<offset)).*t + (t>=offset).*(onset+1);

Ca = (t<=onset).* bathCa + ((t>onset) & (t<offset)).*(bathCa + (((0.12 * I) / (4 * pi * 8e-10 * 1.602e-19 * r))) * erfc(r./(2*(sqrt(8e-10*(terfc-onset))))) / 6.022e23) + (t>=offset).* bathCa


%[Ca] = 0;
%if t<=onset; Ca = bathCa;
%elseif (t>onset) & (t<offset);
%Ca = bathCa + (((0.12 * I) / (4 * pi * 8e-10 * 1.602e-19 * r))) * erfc(r./(2*(sqrt(8e-10*(t-onset))))) / 6.022e23;
%elseif t>=offset;
%Ca(region3) = bathCa;
%end
