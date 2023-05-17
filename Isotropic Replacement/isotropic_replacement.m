function [isotropic_image, calibration_new] = isotropic_replacement(anisotropic_image, calibration_old, calibration_new)
% Correct spatially anisotropic images using isotropic replacement method
% where signal pixels are placed into an isotropic image whose overall
% (real-world) size matches the original image but it has been subdivided
% to make isotropic. Signal pixels are placed into new pixels whose
% (real-world) position most closely matches the position of the original
% pixels. Specify the original pixel size (length-units per pixel) via the
% calibration_old input, and specify the desired new pixel size via the
% (optional) calibration_new inpit. In no new calibration is specified, a
% default size will be used as either the smallest of the original
% calibration values or the largest divided by root-3, whichever value is
% smaller. Input images must be masks in the form of logical matrices. If
% correcting mutliple images, input a cell array of logical matrices. The
% corrected image(s) will be output in the same format as the input (cell
% array of logical matrices VS logical matric). You can specify the desired
% 
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/16/2023

    %% Input Validation

    % check for correct number of inputs
    if nargin < 1
        error("Not enough input arguments.");
    elseif nargin > 3
        error("Too many input arguments.");
    end

    % check shape of calibration input
    if ~isrow(calibration_old)
        error("Calibration must be a row vector.");
    end
    if nargin == 3
        if ~isrow(calibration_old)
            error("Calibration must be a row vector.");
        end
        if ~all(size(calibration_old) == size(calibration_new))
            error("New and old image calibrations must be the same size.");
        end
    end

    % check if first input is correct data types
    is_one_image = false;
    if ~iscell(anisotropic_image)
        is_one_image = true;
        anisotropic_image = {anisotropic_image};
    end
    if ~all(cellfun(@islogical,anisotropic_image))
        error("Images must be a logical matrix or cell array of logical matrices.");
    end

    % check if all images are the same size
    im_sz = cellfun(@size, anisotropic_image, 'UniformOutput',false);
    if numel(unique(cellfun(@numel,im_sz))) > 1
        % check number of dimensions
        error("Images are not all the same size.");
    end
    if ~isrow(unique(cell2mat(im_sz),'rows'))
        % check overall size
        error("Images are not all the same size.");
    end

    % check if image is already isotropic
    if numel(unique(calibration_old)) < 1 
        isotropic_image = anisotropic_image;
        calibration_new = calibration_old;
        warning("Images are already isotropic. Original inputs returned.");
        return;
    end

    %% Correct for Spatial Anisotropy

    % make column vector for consistent formatting
    anisotropic_image = anisotropic_image(:);

    % Calculate new isotropic pixel size if none provided
    if nargin == 2
        min_cal = min(calibration_old);
        max_cal = max(calibration_old);
        root3dim = max_cal / sqrt(3);
        if min_cal <= root3dim
            isotropic_distance = min_cal;
        else
            isotropic_distance = root3dim;
        end
        calibration_new = repelem(isotropic_distance,numel(calibration_old));
    end

    % Get general info
    n_dims = ndims(anisotropic_image{1}); % is the image 2- or 3-D?
    scale_factors = calibration_old ./ isotropic_distance; % how each original calibration will be scaled
    im_sz_old = size(anisotropic_image{1});
    im_sz_new = round(im_sz_old .* scale_factors); % size of corrected image
    n_pixels = numel(anisotropic_image{1}); % calculate the total number of pixels in the original image

    % determine isotropic subscripts for every pixel in the original image
    old_indices = (1:n_pixels)'; % get linear indices of every pixel in the original image
    [y, x, z] = ind2sub(im_sz_old, old_indices); % convert linear indices to subscripts
    new_subscripts = round([y,x,z] .* scale_factors); %convert anisotropic subscripts to isotropic
    new_subscripts = mat2cell(new_subscripts, n_pixels, ones(1,n_dims)); % prep for dimensionality-independent sub2ind conversion
    replacement_indices = sub2ind(im_sz_new, new_subscripts{:}); % convert isotropic subscripts to linear indices
    isotropic_template = false(im_sz_new);

    % correct masks by converting anisotropic indices to isotropic indices
    n_images = numel(anisotropic_image);
    isotropic_image = repelem({false(im_sz_new)}, n_images, 1);
    for i_image = 1:n_images
        isotropic_image{i_image} = replace_pixels(replacement_indices, isotropic_template, anisotropic_image{i_image});
    end

    % return subject mask to matrix if only one was provided
    if is_one_image
        isotropic_image = isotropic_image{1};
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
