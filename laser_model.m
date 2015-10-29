function [pos_resp, neg_resp, total_resp] = laser_model(t,laser_start,laser_end,laser_power)
% This function models the displacement response of a lasered hair bundle.
% The model is based on results from masking experiments, in which it  
% there appears to be a positive and a negative component of the response.
% The positive component is linearly related to power, while the negative
% is exponentially related to power. In this case the fourth power.
% INPUTS:   
%   t           - time vector
%   laser_start - time in seconds at which the laser comes on
%   laser_end   - time in seconds at which laser turns off
%   laser_power - power of the laser
%
% OUTPUTS:  
%   pos_resp    - response given by the positive element
%   neg_resp    - response given by the negative element
%   total_resp  - sum of the positive and negative responses
%
% The sensitivity of the positive response is set by the variable 'sens'
% The sensitivity of the megative response is set by 'm'
%
% --------------------------------- %
%       Julien B. Azimzadeh
%       jazimzadeh@rockefeller.edu
%       September 26, 2015
% --------------------------------- %
%
sens    = 40;                       % 40 nm/mW of laser power
m       = 0.7*(laser_power^4);      % coefficient that confers nonlinearity the negative response's dependence on laser power

% Initialize output vectors
pos_resp = zeros(length(t),1);
neg_resp = zeros(length(t),1);
total_resp = zeros(length(t),1);

for s = 1:length(t)
    if t(s) <= laser_start
        pos_resp(s) = 0;
        neg_resp(s) = 0;
    end
    if t(s) > laser_start && t(s) <= laser_end ;
        pos_resp(s) = sens*laser_power;   % Displacement is in nm; 40 nm/mW
        neg_resp(s) = -m*(t(s)-laser_start).^(1/exp(1));
    end
    if t(s) > laser_end
        pos_resp(s) = 0;
        min = neg_resp(find(laser_end==t));
        neg_resp(s) = min*exp(-15*(t(s)-laser_end));
    end
    total_resp(s) = pos_resp(s) + neg_resp(s);
end