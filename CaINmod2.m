function [CaINmod2] = CaINmod2(time,onset,offset, delta_t,I,r,bathCa)
% onset is the onset of the Ca pulse (in millseconds)
% offset is the offset of the Ca pulse
% duration is the duration of the trace (the [Ca] simply goes back to its
%    bath value 10 ms after the offset
% I is the amplitude of current injected into the pipette, in Amps
% r is the distance from the tip of the pipette to the top of the hair
%   bundle, in meters
%bathCa is the resting concentration of Ca in the bath, given in M

% the concentration of Ca outputed by this function is in mol/liter
concentration = zeros(1,onset/delta_t); %preallocate

pre_ionto = 0 : delta_t : onset;
ionto = onset : delta_t : offset;
post_ionto = offset : delta_t : time;

    for i = 1:1:length(pre_ionto)
        concentration(i) = bathCa;
    end

    for i = 1:1:length(ionto)
%i need i in terms of t
        t = i*delta_t;
        concentration2(i) =  bathCa +  (((0.12 * I) / (4 * pi * 8e-10 * 1.602e-19 * r)) *   erfc(r./(2*(sqrt(8e-10*t)))) / 6.022e23 / 1000);
    end
    
    for i = 1:1:length(post_ionto)
        t=i * delta_t;
        concentration3(i) =  (((0.12 * I) / (4 * pi * 8e-10 * 1.602e-19 * r)) * ( erfc(r./(2*(sqrt(8e-10*t)))) - erfc(r./(2*sqrt((8e-10)*(t-t(1))))) )) /6.022e23/1000;
    end
    subplot(1,3,1)
    plot(concentration)
    subplot(1,3,2)
    plot(concentration2)
    subplot(1,3,3)
    plot(concentration3)
    
end

