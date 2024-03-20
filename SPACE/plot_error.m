function [line_plot, shaded_error] = plot_error(x, y, error_lower, error_upper)
% Function to plot a line with shaded error. Error bounds error_lower and
% error_upper should be absolute error, as opposed to error defined
% relative to the line (y). For example, at point (1,1) on the line, if
% lower and upper error are y+1 and y-1, the corresponding values in
% error_lower and error_upper should be 2 and 0, NOT +1 and -1. 
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% GitHub: https://github.com/andrewsoltisz/SPACE---Spatial-Pattern-Analysis-using-Closest-Events
% Publication: https://doi.org/10.1093/mam/ozae022
% Last Updated: 10/05/2023
%
% Copyright (C) 2024, Andrew Michael Soltisz. All rights reserved.
% This source code is licensed under the BSD-3-Clause License found in the
% LICENSE.txt file in the root directory of this source tree.

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
