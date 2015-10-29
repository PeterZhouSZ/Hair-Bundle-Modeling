function [pos_output] = test_piecewise(t_series,t)
pos_output = zeros(length(t_series),1);

for s = 1:length(t_series)
    if t_series(s) <= t
        pos_output(s) = 0
    end
    if t_series(s) > t & t_series(s) < 8
        pos_output(s) = 1
    end
    if t_series(s) >= 8
        pos_output(s) = 0
    end
end
