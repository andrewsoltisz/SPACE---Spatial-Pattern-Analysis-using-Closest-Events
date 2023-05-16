function [p, h] = ttest2w(x, y, n_x, n_y, alpha)
% Compare means of samples x and y using the 2-sided Student's T-test with
% weighted measurements according to their sample sizes n_x and n_y.
%
% The weighted ttest is computed the same as MATLAB's built-in function
% 'ttest2', but the weighted mean and variance are used in place of their
% unweighted counterparts. Weighted values are defined on the webpage:
% https://support.sas.com/documentation/onlinedoc/stat/131/ttest.pdf
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/11/2023

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
n_x = n_x(:);
n_y = n_y(:);

% degrees of freedom (df)
sz_x = numel(x);
sz_y = numel(y);
dfe = sz_x + sz_y - 2;

% weighted mean
mu_w = @(a, n_a) sum(a .* n_a) / sum(n_a);
mu_w_x = mu_w(x, n_x);
mu_w_y = mu_w(y, n_y);

% weighted variance
var_w = @(a, w, n, mu_w) sum(n.*((a - mu_w).^2)) / (n - 1);
var_w_x = var_w(x, n_x, sz_x, mu_w_x);
var_w_y = var_w(y, n_y, sz_y, mu_w_y);
difference = mu_w_x - mu_w_y;

% Calculate p-value
sPooled = sqrt(((sz_x - 1) .* var_w_x + (sz_y - 1) .* var_w_y) ./ dfe);
se = sPooled .* sqrt((1 ./ sz_x) + (1 ./ sz_y));
ratio = difference ./ se;
p = 2 * tcdf(-abs(ratio), dfe);
h = cast(p <= alpha, 'like', p);
h(isnan(p)) = NaN; 

end
