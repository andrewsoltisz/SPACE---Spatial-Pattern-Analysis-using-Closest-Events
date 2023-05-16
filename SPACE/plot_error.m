function [line_plot, shaded_error] = plot_error(x, y, error_lower, error_upper)
% Function to plot a line with shaded error. If error is symetric about the
% line, then specify y-coordinates of the absolute error (as opposed to
% error relative to the line) as a single third input. If error is
% asymmetric about the line with distinct upper and lower regimes, then
% specify the y-coordinates of each regime as a third (lower) and fourth
% (upper) inputs.
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/11/2023

    % Input validation
    if nargin < 4
        error("Not enough input arguments.");
    elseif nargin > 4
        error("Too many input arguments");
    end
    if numel(unique(arrayfun(@numel, x, y, error_lower, error_upper))) > 1
        error("Input sizes do not match.");
    end

    % Line with shaded error 
    hold on;
    shaded_error = patch('XData', [x;flip(x)], 'YData', [error_lower;flip(error_upper)], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    line_plot = plot(x, y);
    hold off;

end
