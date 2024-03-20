function [p, h] = ttest2w(x, y, w_x, w_y, alpha)
% Compare means of samples x and y using the 2-sided Student's T-test with
% weighted measurements according to their sample sizes w_x and w_y.
%
% The weighted ttest is computed the same as MATLAB's built-in function
% 'ttest2', but the weighted mean and variance are used in place of their
% unweighted counterparts.
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

if nargin < 5
    alpha = 0.05; % default value
else
    if (alpha > 1 || alpha < 0) || ~isscalar(alpha)
        error("Alpha must be a scalar number between 0 and 1.");
    end
end

% convert to column vector
x = x(:);
y = y(:);
w_x = w_x(:);
w_y = w_y(:);

% degrees of freedom (df)
n_x = numel(x);
n_y = numel(y);
dfe = n_x + n_y - 2;

% weighted variance
[var_w_x, mu_w_x] = var(x, w_x);
[var_w_y, mu_w_y] = var(y, w_y);
mu_difference = mu_w_x - mu_w_y;

% Calculate p-value (based on code found in built-in ttest2)
sPooled = sqrt(((n_x - 1) .* var_w_x + (n_y - 1) .* var_w_y) ./ dfe);
se = sPooled .* sqrt((1 ./ n_x) + (1 ./ n_y));
ratio = mu_difference ./ se;
p = 2 * tcdf(-abs(ratio), dfe);
h = cast(p <= alpha, 'like', p);
h(isnan(p)) = NaN; 

end
