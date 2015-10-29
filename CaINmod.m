function [CaINmod] = CaINmod(time,onset,offset, delta_t,I,r,bathCa)
% onset is the onset of the Ca pulse (in millseconds)
% offset is the offset of the Ca pulse
% duration is the duration of the trace (the [Ca] simply goes back to its
%    bath value 10 ms after the offset
% I is the amplitude of current injected into the pipette, in Amps
% r is the distance from the tip of the pipette to the top of the hair
%   bundle, in meters
%bathCa is the resting concentration of Ca in the bath, given in M

% the concentration of Ca outputed by this function is in mol/liter

tspan = 0:delta_t:time;

for i = 1:1:length(tspan)
    t = tspan(i);         % time vector of correct duration
    C(i) =    (((0.12 * I) / (4 * pi * 8e-10 * 1.602e-19 * r)) *   erfc(r./(2*(sqrt(8e-10*t)))) / 6.022e23 / 1000);
    Coff(i) = (((0.12 * I) / (4 * pi * 8e-10 * 1.602e-19 * r)) * ( erfc(r./(2*(sqrt(8e-10*t)))) - erfc(r./(2*sqrt((8e-10)*(t-tspan(1))))) )) /6.022e23/1000;

    if (i < onset)
        concentration(i) = bathCa;
    end

    if (i > onset & i < offset)
        %C_down = downsample(C,100);
        concentration(i) = bathCa + C(i-onset/delta_t);
    end
    
    if (i >= offset )   
        %Coff_down = downsample(Coff,100);
        concentration(i) = bathCa + Coff(i-offset/delta_t);
    end
    
end
plot(tspan,concentration)
