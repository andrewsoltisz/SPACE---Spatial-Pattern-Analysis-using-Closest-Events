function [isotropic_mask, isotropic_ROI, calibration_new] = isotropic_replacement(anisotropic_mask, anisotropic_ROI, calibration_old, calibration_new)
% Resample spatially anisotropic images using isotropic replacement method
% where signal (TRUE) pixels are placed into an isotropic image whose overall
% (real-world) size matches the original input image but has been subdivided
% to be isotropic. Signal pixels are then placed into new pixels whose
% (real-world) position most closely matches the position of the original
% pixels. 
% 
% -------------------------------------------------------------------------
% 
% INPUTS: 
% 
% 1. anisotropic_mask - Binary image (mask) identifying the positions of
%                       events or signals in the original (anisotropic)
%                       image. Here, elements with a value of 1 (or TRUE)
%                       indicate pixels containing signal and elements with
%                       a value of 0 (or FALSE) indicate background pixels.
%                       When resampling only one mask, format as a logical
%                       matrix of any dimension with the same size and
%                       shape as 'anisotropic_ROI'. To resample multiple
%                       masks with one function call to
%                       isotropic_replacement, format as a cell array of
%                       logical matrices, where all matrices have the same
%                       size and shape.
%
% 2. anisotropic_ROI - Binary image (mask) identifying a region-of-interest
%                      (ROI) within 'anisotropic_mask' where the resampling
%                      will be focused. This input is mandatory, but if you
%                      don't know what it should be, just input a logical
%                      matrix of all 1's. Here, elements with a value of 1
%                      (or TRUE) indicate pixels containing signal and
%                      elements with a value of 0 (or FALSE) indicate
%                      background pixels. Formatted as a logical matrix of
%                      any dimension with the same size and shape as the
%                      matrix component of 'anisotropic_mask'.
%
% 3. calibration_old - Specifies the real-world dimensions (side-lengths)
%                      of pixels composing 'anisotropic_mask' and
%                      'anisotropic_ROI'. Formatted as a 1-by-d numerical
%                      row vector where each column specifies the pixel
%                      side-length along each dimension in the order x, y,
%                      z, etc. All masks should have the same original
%                      calibration, even when resampling multiple images
%                      with one function call to isotropic_replacement.
%
% 4. calibration_new - OPTIONAL input speicifying the final cuboidal
%                      side-length of the pixels which the output masks
%                      will be resample to match. Formatted as a positive
%                      numeric scalar that is smaller than or equal the
%                      smallest value in 'calibration_old'. If this input
%                      is not provided, a default value will be used that
%                      is either the smallest value from 'calibration_old'
%                      or the largest value from 'calibration_old' divided
%                      by root-3, whichever value is smallest. This latter
%                      value is chosen to ensure that any resulting errors
%                      on pair-wise distances between pixels is less than
%                      the lowest-resolution dimension from the original
%                      image.
% 
% -------------------------------------------------------------------------
% 
% OUTPUTS: 
% 
% 1. isotropic_mask - Binary image that is a resampling of
%                     'anisotropic_mask' such that its pixels are isotropic
%                     or cuboidal. Formatted the same as
%                     'anisotropic_mask'.
%
% 2. isotropic_ROI - Binary image that is a resampling of
%                     'anisotropic_ROI' such that its pixels are isotropic
%                     or cuboidal. Formatted the same as 'anisotropic_ROI'.
%
% 3. calibration_new - Specifies the real-world size (side-length)
%                      of the cuboidal pixels composing 'isotropic_mask'
%                      and 'calibration_new'. Foratted as a positive
%                      numeric scalar. 
% -------------------------------------------------------------------------
% 
% USE CASES: 
% 
% [...] = isotropic_replacement(anisotropic_mask, anisotropic_ROI, calibration_old)
% perform Isotropic Replacement on 'anisotropic_mask' and
% 'anisotropic_ROI'. The new pixel size 'calibration_new' will be
% automatically calculated as either the smallest value from
% 'calibration_old' or the largest value from 'calibration_old' divided by
% root-3, whichever value is smallest.
%
% [...] isotropic_replacement(anisotropic_mask, anisotropic_ROI, calibration_old, calibration_new)
% % perform Isotropic Replacement on 'anisotropic_mask' and
% 'anisotropic_ROI'. The new pixel size 'calibration_new' will be defined
% by user input. 
%
% -------------------------------------------------------------------------
% 
% AUTHORSHIP: 
% 
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% GitHub: https://github.com/andrewsoltisz/SPACE---Spatial-Pattern-Analysis-using-Closest-Events
% Publication: https://doi.org/10.1101/2023.05.17.541131
% Last Updated: 10/05/2023
%
% Copyright 2023, Andrew Michael Soltisz. All rights reserved.
% This source code is licensed under the BSD-3-Clause License found in the
% LICENSE.txt file in the root directory of this source tree.
% 
% -------------------------------------------------------------------------

    %% Input Validation

    % check for correct number of inputs
    if nargin < 1
        error("Not enough input arguments.");
    elseif nargin > 4
        error("Too many input arguments.");
    end

    % check shape of calibration input
    if ~isrow(calibration_old)
        error("Calibration must be a row vector.");
    end
    if nargin == 4
        if ~isrow(calibration_old)
            error("Calibration must be a row vector.");
        end
        % make sure new calibration is isotropic
        if numel(calibration_new) > 1 || ~isnumeric(calibration_new) || calibration_new <= 0
            error("New calibration must be a positive numeric scalar.");
        end
        % make sure new calibration is no more than old calibration
        if any(calibration_new > calibration_old)
            error("New calibration must be <= old calibration.");
        end
    end

    % check if first input is correct data types
    just_one_image = false;
    if ~iscell(anisotropic_mask)
        just_one_image = true;
        anisotropic_mask = {anisotropic_mask};
    end
    if ~all(cellfun(@islogical,anisotropic_mask))
        error("Masks must be a logical matrix or cell array of logical matrices.");
    end
    if ~islogical(anisotropic_ROI)
        error("ROI mask must be a logical matrix.");
    end


    % check if all images are the same size
    im_sz = cellfun(@size, anisotropic_mask, 'UniformOutput',false);
    if numel(unique(cellfun(@numel,im_sz))) > 1
        % check number of dimensions
        error("Masks are not all the same size.");
    end
    if ~isrow(unique(cell2mat(im_sz),'rows'))
        % check overall size
        error("Masks are not all the same size.");
    end
    im_sz_old = size(anisotropic_mask{1});
    if ~all(im_sz_old == size(anisotropic_ROI))
        error("ROI mask size does not match other masks.");
    end

    % check if image is already isotropic
    if numel(unique(calibration_old)) < 1 
        isotropic_mask = anisotropic_mask;
        calibration_new = calibration_old;
        warning("Images are already isotropic. Original inputs returned.");
        return;
    end

    %% Correct for Spatial Anisotropy

    % make column vector for consistent formatting
    anisotropic_mask = anisotropic_mask(:);

    % Calculate new isotropic pixel size if none provided
    if nargin == 3
        [min_cal, min_idx] = min(calibration_old);
        [~, max_idx] = max(calibration_old);
        sqrt3_cal = calibration_old / sqrt(3);
        if min_cal <= sqrt3_cal(max_idx) % no error introduced to min_cal dimensions, all other error is sub-resolution
            isotropic_distance = min_cal;
        else 
            isotropic_distance = sqrt3_cal(min_idx);
        end
        calibration_new = isotropic_distance;
    elseif nargin == 4
        isotropic_distance = calibration_new;
    end

    % Get general info
    n_dims = ndims(anisotropic_mask{1}); % is the image 2- or 3-D?
    scale_factors = calibration_old ./ isotropic_distance; % how each original calibration will be scaled
    im_sz_new = round(im_sz_old .* scale_factors); % size of corrected image
    n_pixels = numel(anisotropic_mask{1}); % calculate the total number of pixels in the original image

    % determine isotropic subscripts for every pixel in the original image
    old_indices = (1:n_pixels)'; % get linear indices of every pixel in the original image
    [y, x, z] = ind2sub(im_sz_old, old_indices); % convert linear indices to subscripts
    new_subscripts = round([y,x,z] .* scale_factors); %convert anisotropic subscripts to isotropic
    new_subscripts = mat2cell(new_subscripts, n_pixels, ones(1,n_dims)); % prep for dimensionality-independent sub2ind conversion
    replacement_indices = sub2ind(im_sz_new, new_subscripts{:}); % convert isotropic subscripts to linear indices
    isotropic_template = false(im_sz_new);

    % correct masks by converting anisotropic indices to isotropic indices
    n_images = numel(anisotropic_mask);
    isotropic_mask = repelem({false(im_sz_new)}, n_images, 1);
    for i_image = 1:n_images
        isotropic_mask{i_image} = replace_pixels(replacement_indices, isotropic_template, anisotropic_mask{i_image});
    end
    isotropic_ROI = replace_pixels(replacement_indices, isotropic_template, anisotropic_ROI);

    % return subject mask to matrix data-type if only one was input
    if just_one_image
        isotropic_mask = isotropic_mask{1};
    end
    
end

%% Helper Functions

function [isotropic_mask, isotropic_indices] = replace_pixels(replacement_indices, isotropic_mask, anisotropic_mask)
% Given isotropic indices for every pixel in the original image
% (replacement_indices), place the true values from the anisotropic_mask
% into the pixels of isotropic_mask closest to their origina (real-world)
% positions.
% 
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/16/2023

    isotropic_indices = replacement_indices(anisotropic_mask(:));
    isotropic_mask(isotropic_indices) = true;

end
