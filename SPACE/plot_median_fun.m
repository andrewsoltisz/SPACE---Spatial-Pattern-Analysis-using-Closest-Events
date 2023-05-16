function [median_line, quantile_envelope] = plot_median_fun(x, y_med, y_lower, y_upper, c, ls)
% Function to plot a median empirical function (x, y_med) with shaded quantile
% envelopes bounded by the lower (x, y_lower) and upper (x, y_upper)
% quantiles. Optional arguments can be added to specify the line/envelope
% color (c) and linestyle (ls). 
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/11/2023

    hold on;
    quantile_envelope = patch('XData', [x;flip(x)], 'YData', [y_lower;flip(y_upper)], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    median_line = plot(x, y_med);

    if nargin > 4
        quantile_envelope.FaceColor = c;
        median_line.Color = c;
    end
    if nargin > 5
        median_line.LineStyle = ls;
    end
    hold off;

end
