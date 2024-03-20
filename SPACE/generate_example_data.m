% Generate synthetic data example scripts that illustrate how to use SPACE.
% 3 groups of images are generated, each with distinct X<->Y spatial
% relationships: dispersed (dis), independent (ind), and aggregated (agg).
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% GitHub: https://github.com/andrewsoltisz/SPACE---Spatial-Pattern-Analysis-using-Closest-Events
% Publication: https://doi.org/10.1101/2023.05.17.541131
% Last Updated: 10/05/2023
%
% Copyright (C) 2024, Andrew Michael Soltisz. All rights reserved.
% This source code is licensed under the BSD-3-Clause License found in the
% LICENSE.txt file in the root directory of this source tree.

%% Prepare environment

close all; 
clear all; 
clc;

%% Define image parameters

n_images = 50; % number of images per group
n_total = n_images * 3; % total number of images across groups
mask_size = [100, 100]; 
roi_mask = true(mask_size); % define region-of-interest (ROI) mask as full image 
ran_n = prod(mask_size); % number of pixels composing the ROI and therefore the number of events used to generate the random distribution 
% create some variability for more interesting data 
X_conc = abs(normrnd(0.02, 0.001, n_images, 1)); % event concentration of pattern-X 
Y_conc = abs(normrnd(0.02, 0.001, n_images, 1)); % event concentration of pattern-Y 
% Y-->X spatial proximity (S)
S_dis = normrnd(-0.55, 0.01, n_images, 1); % spatially dispersed 
S_ind = normrnd(0, 0.01, n_images, 1); % spatially independent 
S_agg = normrnd(0.8, 0.01, n_images, 1); % spatially aggregated 
% plotting colors
X_color = 'r';
Y_color = 'g';

%% Preallocate images
X_mask_list_dis = repelem({false(mask_size)}, n_images, 1);
Y_mask_list_dis = X_mask_list_dis;
X_mask_list_ind = X_mask_list_dis;
Y_mask_list_ind = X_mask_list_dis;
X_mask_list_agg = X_mask_list_dis;
Y_mask_list_agg = X_mask_list_dis;
X_n_dis = nan(n_images, 1);
Y_n_dis = X_n_dis;
X_n_ind = X_n_dis;
Y_n_ind = X_n_dis;
X_n_agg = X_n_dis;
Y_n_agg = X_n_dis;
overlay_dis = repelem({zeros([mask_size, 3])}, n_images, 1); 
overlay_ind = overlay_dis;
overlay_agg = overlay_dis;
ROI_mask_list_dis = repelem({roi_mask},n_images);
ROI_mask_list_ind = ROI_mask_list_dis;
ROI_mask_list_agg = ROI_mask_list_dis;

%% Generate Images

wb = waitbar(0, sprintf("Generating image: %i/%i", 0, n_total));
for i_image = 1:n_images

    % create image masks
    [Y_mask_list_dis{i_image}, X_mask_list_dis{i_image}] = gen_synthetic_masks(mask_size, X_conc(i_image), Y_conc(i_image), S_dis(i_image));
    [Y_mask_list_ind{i_image}, X_mask_list_ind{i_image}] = gen_synthetic_masks(mask_size, X_conc(i_image), Y_conc(i_image), S_ind(i_image));
    [Y_mask_list_agg{i_image}, X_mask_list_agg{i_image}] = gen_synthetic_masks(mask_size, X_conc(i_image), Y_conc(i_image), S_agg(i_image));

    % update waitbar
    waitbar(i_image/n_images, wb, sprintf("Generating image: %i/%i", i_image*3, n_total));
end
close(wb);

%% Save Data For example_single

example_data = 'Example Data';
if ~isfolder(example_data)
    mkdir(example_data);
end
example_data = [example_data,filesep];

imwrite(X_mask_list_agg{1}, [example_data, 'X_mask.tif']);
imwrite(Y_mask_list_agg{1}, [example_data, 'Y_mask.tif']);
imwrite(ROI_mask_list_agg{1}, [example_data, 'ROI_mask.tif']);

%% Save Data For example_batch

% create clearer variable names for example script
X_mask_list = X_mask_list_dis;
Y_mask_list = Y_mask_list_dis;
ROI_mask_list = ROI_mask_list_dis;

save([example_data,'X_mask_list.mat'], 'X_mask_list');
save([example_data,'Y_mask_list.mat'], 'Y_mask_list');
save([example_data,'ROI_mask_list.mat'], 'ROI_mask_list');

%% Save Data For example_compare_groups

% Dispersed
save([example_data,'X_mask_list_dis.mat'], 'X_mask_list_dis');
save([example_data,'Y_mask_list_dis.mat'], 'Y_mask_list_dis');
save([example_data,'ROI_mask_list_dis.mat'], 'ROI_mask_list_dis');

% Independent
save([example_data,'X_mask_list_ind.mat'], 'X_mask_list_ind');
save([example_data,'Y_mask_list_ind.mat'], 'Y_mask_list_ind');
save([example_data,'ROI_mask_list_ind.mat'], 'ROI_mask_list_ind');

% Aggregated
save([example_data,'X_mask_list_agg.mat'], 'X_mask_list_agg');
save([example_data,'Y_mask_list_agg.mat'], 'Y_mask_list_agg');
save([example_data,'ROI_mask_list_agg.mat'], 'ROI_mask_list_agg');

%% Helper Functions

function [Y_im, X_im, Y_n, X_n] = gen_synthetic_masks(im_sz, X_conc, Y_conc, S, X_im)
% Generate sythetic images depicting 2 spatial patterns, X and Y, such that
% pattern Y's spatial proximity is tuned by a single input parameter (S).
% Given an specified image size (im_sz) and concentrations of pattern X
% (X_conc) and Y (Y_conc), X-events are assigned to independent and ranom
% pixels within the image afterwhich Y-events are assigned to specific
% pixels whose distance-from-X are chosen using a custom, inverse-CDF-based
% method based on the value of S whose viable range is [-1,1]. The more
% negative S is, the more dispersed Y will be relative to X. The more
% positive S is, the more aggregated Y will be relative to X. And the
% closer to 1 S is, the more independent Y will be relative to X. When
% S==0, Y-events will be assigned to random pixels completely independent
% of X-events. When S==1, all Y-events will be assigned to pixels occupied
% by X-events, ensureing 100% Y-X overlap. When S==-1, Y-events will be
% assigned to the pixel(s) that are furthest from X-events. Intermediate
% values of S will produce intermediate spatial relationships. 
%
% Warning: Y-events are assigned to pixels using a memory-inefficient
% vectorized algorithm, so this algorithm will likely only support the
% creation of images with fewer than 1e6 pixels before RAM overload, though
% this will be system-dependent. 
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/11/2023

    im_n_pixels = prod(im_sz); % calculate total number of pixels in the image
    
    %% Define X-Events and create X-Mask
    
    if nargin < 5 % make a X-image if one isn't provided
        X_n_points = round(X_conc * im_n_pixels); % calculate number of X-events based on X's concentration
        X_idx = randperm(im_n_pixels, X_n_points); % select random linear indices of pixels for X-events to be assigned
        X_im = false(im_sz); % create blank X-mask
        X_im(X_idx) = true; % place X-events into the X-mask
        X_n = sum(X_im, 'all'); % count number of X-pixels
    else
        if any(im_sz ~= size(X_im))
            error("Image size does not match input X-image size.");
        end
    end

    %% Generate Distribution for Assigning Distances to Y-Events
    
    X_dt = bwdist(X_im); % generate distance transformation of X-image
    [X_pdf_x, Y_pdf_y] = epdf(X_dt(:)); % generate distance PDF
    
    filter_y = pdf_filter(S, X_pdf_x); % generate filter: half-normal
    X_pdf_y_filt = filter_y .* Y_pdf_y'; % apply filter to PDF y-coords
    X_cdf_y_filt = cumsum(X_pdf_y_filt); % generate CDF y-values
    X_cdf_y_filt = X_cdf_y_filt ./ X_cdf_y_filt(end);

    %% Define Y-Events and Create Y-Mask
    
    Y_n_points = round(Y_conc * im_n_pixels); % Calculate number of Y-events based on its concentration

    % Interp1 (used below) requires all sample values to be unqiue.
    % Adjacent uniques are used to avoid skipping over the curved
    % transition between 0 and 1 for CDFs with sudden changes (when S is
    % close to +/- 1). This transition is on the right ('last') when S<0
    % and on the left ('first') when S>0.
    if S < 0
        [X_cdf_y_uniq, idxUnique] = unique(X_cdf_y_filt', 'last');
    else
        [X_cdf_y_uniq, idxUnique] = unique(X_cdf_y_filt', 'first');
    end

    % Inverse (CDF) transformation method to assign a distance to each Y-event
    X_cdf_x_iF = X_pdf_x(idxUnique);
    Y_cdf_y = (X_cdf_y_uniq(end) - X_cdf_y_uniq(1)) .* rand([1,Y_n_points]) +  X_cdf_y_uniq(1); % give each Y-event a random CDF y-value within the X_cdf_y range
    if numel(X_cdf_y_uniq) == 1
        Y_cdf_x = repelem(X_cdf_x_iF, Y_n_points);
    else
        Y_cdf_x = interp1(X_cdf_y_uniq, X_cdf_x_iF, Y_cdf_y); % determine x-value for each Y-event CDF y-value
    end

    % Find random pixels with nearest distance to assigned Y-distances
    [X_dt_sorted_vals, X_dt_sorted_idx] = sort(X_dt(:)); % sort distances for vectorized computations later on
    [distUnique, distIdx] = unique(X_dt_sorted_vals); % determine all unique distance measurements and their indices
    [~, Y_distIdx] = min(abs(Y_cdf_x - distUnique), [], 1); % find index of closest distance to each assigned Y-distance
    idxLow = distIdx(Y_distIdx); % find lower bound of pixel indices that Y-events can be assigned
    edgeCase = Y_distIdx == numel(distIdx); % these Y-points are in last idx bin and getting an idx range for them (3 lines below) will throw out-of-range error. Handle separately.
    idxHigh = zeros(Y_n_points, 1); % preallocate
    idxHigh(edgeCase) = im_n_pixels; % points in last idx bin have upper idx bound equal to total number of pixels
    idxHigh(~edgeCase) = distIdx(Y_distIdx(~edgeCase) + 1) - 1; % find upper bound of pixel indices that Y-events (with multiple options) can be assigned
    Y_idx = X_dt_sorted_idx(round((idxHigh - idxLow) .* rand([Y_n_points,1]) + idxLow)); % randomly assign each Y-event to a pixel within its idx range
   
    % Populate Y-Mask
    Y_im = false(im_sz); % create blank Y-mask
    Y_im(Y_idx) = true; % place Y-events into Y-mask
    Y_n = sum(Y_im, 'all'); % count number of Y-pixels

end

function filter_y = pdf_filter(S, x_range)
% Define function to convert S to the standard deviation (sigma) of the
% half-normal function used to filter pattern Y's inverse-CDF
% distribution. 
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/11/2023

    res = 100; % resolution (res) of transfer function
    transfer_y = logspace(-1.585, 0.60206, res); 
    transfer_x = linspace(1,0,res);
    S_abs = abs(S);
    sigma = interp1(transfer_x, transfer_y, S_abs, 'pchip'); % standard deviation of filtering half-normal distribution

    % Preallocate
    n_x = numel(x_range);
    filter_y = zeros(1, n_x);
    filter_x = linspace(0, 1, n_x);
    
    % Generate filter based on half-normal distribution
    if S_abs == 1 % Degenerate cases: S==+1 or S==-1
        filter_y(1) = 1;
    elseif S == 0 % complete spatial randomness
        filter_y(:) = 1;
    else % intermediate spatial relationships
        filter_y = pdf(makedist('HalfNormal', 'mu', 0, 'sigma', sigma), filter_x);
    end

    filter_y = filter_y ./ max(filter_y); % normalize to max
    
    % The half-normal filter for negative S are mirror of their positive S
    % filter
    if S < 0 
        filter_y = flip(filter_y); 
    end

end

function [pdf_x, pdf_y] = epdf(distances, weight)
% Calculate the empirical probability density function (ePDF) of a set of
% distance measurements. Optionally, events can be weighted, if for example
% each event occupies a volume of space.
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/11/2023

    [pdf_y,pdf_x] = groupcounts(distances);

    % make sure PDF is defined for full distance range [0, maxDistance]
    if pdf_x(1) ~= 0
        pdf_x = [0; pdf_x];
        pdf_y = [0; pdf_y];
    end

    % scale event counts if their weight is specified
    if nargin == 3 
        pdf_y = pdf_y * weight;
    end

end
