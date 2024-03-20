function [im] = gen_overlay(X_mask,Y_mask)
% Create an RGB overlay image of 2 masks by placing the first mask in the red
% channel and the second in the green channel.
%
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% GitHub: https://github.com/andrewsoltisz/SPACE---Spatial-Pattern-Analysis-using-Closest-Events
% Publication: https://doi.org/10.1093/mam/ozae022
% Last Updated: 10/05/2023
%
% Copyright (C) 2024, Andrew Michael Soltisz. All rights reserved.
% This source code is licensed under the BSD-3-Clause License found in the
% LICENSE.txt file in the root directory of this source tree.

    X_sz = size(X_mask);
    Y_sz = size(Y_mask);

    if ~all(X_sz == Y_sz)
        error("Image mask sizes must match.");
    end

    im = zeros([X_sz,3]);
    im(:,:,1) = X_mask;
    im(:,:,2) = Y_mask;

end
