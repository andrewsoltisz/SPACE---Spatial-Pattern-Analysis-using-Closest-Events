function [im] = gen_overlay(X_mask,Y_mask)
% Create an RGB overlay image of 2 masks by placing the first mask in the red
% channel and the second in the green channel.
%
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/09/2023

    X_sz = size(X_mask);
    Y_sz = size(Y_mask);

    if ~all(X_sz == Y_sz)
        error("Image mask sizes must match.");
    end

    im = zeros([X_sz,3]);
    im(:,:,1) = X_mask;
    im(:,:,2) = Y_mask;

end