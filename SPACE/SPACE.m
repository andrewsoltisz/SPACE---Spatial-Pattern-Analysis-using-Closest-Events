function [Results_X, Results_Y] = SPACE(X_mask, Y_mask, ROI_mask, pixel_size)
% Perform 'Spatial Pattern Analysis using Closest Events' (SPACE), a
% discrete, nearest neighbor-based point patten analysis which
% characterizes the spatial association between imaged patterns X and Y
% whose events are marked by the 'true' pixels in their corresponding
% binary image masks 'X_mask' and 'Y_mask'.
% 
% -------------------------------------------------------------------------
%
% INPUTS:
%
% 1. X_mask - Binary image (mask) identifying the first signal or type of
%             objects to localize relative to the second signal represented
%             by 'Y_mask'. Here, elements with a value of 1 (or TRUE)
%             indicate pixels containing signal and elements with a value
%             of 0 (or FALSE) indicate background pixels. Formatted as a
%             logical matrix of any dimensionality with the same size and
%             shape as 'Y_mask'.
%
% 2. Y_mask - Binary image (mask) identifying the second signal or type of
%             objects to localize relative to the first signal represented
%             by 'X_mask'. Here, elements with a value of 1 (or TRUE)
%             indicate pixels containing signal and elements with a value
%             of 0 (or FALSE) indicate background pixels. Formatted as a
%             logical matrix of any dimensionality with the same size and
%             shape as 'X_mask'.
%
% 3. ROI_mask - OPTIONAL input which is a binary image (mask) identifying a
%               region-of-interest (ROI) within the other masks where the
%               analysis will be focused. If this is not provided or is
%               assigned the value of an empty matrix, then the intire mask
%               will be analyzed. Here, elements with a value of 1 (or
%               TRUE) indicate the pixels to include in the analysis and
%               elements with a value of 0 (or FALSE) indicate pixels to
%               omit from the analysis. Formatted as a logical matrix of
%               any dimensionality with the same size and shape as 'X_mask'
%               and 'Y_mask'. This input is MANDATORY if 'X_mask' and
%               'Y_mask' were first processed using the Isotropic
%               Replacement method, otherwise this analyses will produce
%               100% garbage results.
%
% 4. pixel_size - OPTIONAL input which specifies the real-world size
%                 (side-length) of pixels composing 'X_mask' and 'Y_mask'.
%                 Formatted as a positive numeric scalar. If this input is
%                 not provided, then distance-derived SPACE outputs will
%                 have units of pixel-count. This input is used to scale
%                 SPACE results so they have the same-real world units as
%                 the input masks. For the SPACE analysis to produce
%                 results which are not 100% garbage, the real-world shape
%                 of mask pixels must be cuboidal, meaning every side of
%                 the pixel is the same length, and this length is the
%                 single value that needs to be assigned to this fouth
%                 input. If your pixels are in fact not cuboidal (e.g. the
%                 x/y-resolution of your masks is different than their
%                 z-resolution), fear not! You just have to resample the
%                 masks so that they are isotropic and have cuboidal
%                 pixels. This can easily be done using the Isotropic
%                 Replacement method found at the GitHub link below. The
%                 'isotropic_ROI mask' (2nd) output from isotropic
%                 replacement function must then be used as the 3rd input
%                 to this function ('ROI_mask') for the SPACE analysis to
%                 be performed correctly.
%
% -------------------------------------------------------------------------
%
% OUTPUTS: 
%
% 1. Results_X - Table containing results for the SPACE analysis of signal
%                identified by X_mask relative to Y_mask. In other words,
%                the analysis is performed using the nearest neighbor
%                distances for all X_mask pixels. Formatted as a table with
%                1 row and 12 columns corresponding to each SPACE output.
%                See the linked GitHub page below for info on what each
%                column represents. If at least one of the X or Y masks
%                contains no signal, then analysis results will be
%                undefined and most table fields will be set to 'NaN'.
%                Results from both outputs are needed to fully describe the
%                spatial assocation between the signals iin X_mask and
%                Y_mask.
%
% 2. Results_Y - Table containing results for the SPACE analysis of signal
%                identified by Y_mask relative to X_mask. In other words,
%                the analysis is performed using the nearest neighbor
%                distances for all Y_mask pixels. Formatted as a table with
%                1 row and 12 columns corresponding to each SPACE output.
%                See the linked GitHub page below for info on what each
%                column represents. If at least one of the X or Y masks
%                contains no signal, then analysis results will be
%                undefined and most table fields will be set to 'NaN'.
%                Results from both outputs are needed to fully describe the
%                spatial assocation between the signals iin X_mask and
%                Y_mask.
%
% -------------------------------------------------------------------------
%
% USES CASES:
% 
% [...] = SPACE(X_mask, Y_mask) 
% Perform analysis without specifiying an ROI or pixel size. Here, the ROI
% is treated as the entire image and the pixel is given a size (side
% length) of 1.
%
% [...] = SPACE(X_mask, Y_mask, ROI_mask) 
% Perform analysis of subregions specified by ROI_mask. Pixels are given a
% size of 1.
%
% [...] = SPACE(X_mask, Y_mask, [], pixel_size) 
% Perform analysis on entire image but pixels are given sizes according to
% pixel_size.
%
% [...] = SPACE(X_mask, Y_mask, ROI_mask, pixel_size) 
% Perform subregion analysis with user defined pixel size.
%
% [Results_X] = SPACE(...) 
% Only return results for analysis of X_mask relative to Y_mask. The other
% analysis (Y_mask relative to X_mask) will not be performed to reduce
% computational load.
%
% [Results_X, Results_Y] = SPACE(...)
% Return results for both analyses - X_mask relative to Y_mask; Y_mask
% relative to X_mask. 
% -------------------------------------------------------------------------
%
% AUTHORSHIP:
%
% Title: Spatial Pattern Analysis using Closest Events (SPACE)
% Citation: Andrew M Soltisz, Peter F Craigmile, Rengasayee Veeraraghavan.
%           Spatial Pattern Analysis using Closest Events (SPACE)—A Nearest Neighbor
%           Point Pattern Analysis Framework for Assessing Spatial Relationships from
%           Digital Images. Microscopy and Microanalysis, 2024, ozae022,
%           https://doi.org/10.1093/mam/ozae022.
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% GitHub: https://github.com/andrewsoltisz/SPACE---Spatial-Pattern-Analysis-using-Closest-Events
% Last Updated: 03/20/2024
%
% Copyright 2024, Andrew Michael Soltisz. All rights reserved.
% This source code is licensed under the BSD-3-Clause License found in the
% LICENSE.txt file in the root directory of this source tree.
%
% -------------------------------------------------------------------------


    % Input validation **************

    % Check for appropriate number of inputs
    if nargin < 2
        error("Not enough input arguments.");
    end

    % Check for appropriate mask data type
    if ~islogical(X_mask)
        error("First input mask is an invalid data type.");
    end
    if ~islogical(Y_mask)
        error("Second input mask is an invalid data type.");
    end

    X_size = size(X_mask);
    Y_size = size(Y_mask);

    % Define default ROI
    if nargin < 3 % no ROI provided, create default
        ROI_mask = true(X_size);
    elseif nargin > 2
        if isempty(ROI_mask) % ROI is empty, create default
            ROI_mask = true(X_size);
        end
    end
    ROI_size = size(ROI_mask);

    % Make sure masks are all the same size
    if ~all(X_size == Y_size) || ~all(X_size == ROI_size)
        error("Mask sizes do not match.");
    end

    % Make sure all masks are logical matrices
    if ~islogical(X_mask) || ~islogical(Y_mask) || ~islogical(ROI_mask)
        error("Masks must be logical matrices.");
    end

    % Make sure pixel length 
    if nargin == 4
        if ~(isscalar(pixel_size) && isnumeric(pixel_size))
            error("Pixel size must be scalar.");
        end
    else
        pixel_size = 1; % default value
    end


    % Compute distribution functions  **************

    % check if both masks contain signal
    X_mask_hasSignal = sum(X_mask,'all') > 0;
    Y_mask_hasSignal = sum(Y_mask,'all') > 0;
    bothMasksHaveSignal = X_mask_hasSignal & Y_mask_hasSignal;
   
    % X-->Y analysis
    if bothMasksHaveSignal
        Results_X = SPACE_normalOutput(X_mask, Y_mask, ROI_mask, pixel_size);
    else
        Results_X = SPACE_nullOutput(X_mask, ROI_mask);
    end

    % Y-->X analysis
    if nargout == 2 % skip analysis if user doesn't ask for it
        if bothMasksHaveSignal
            Results_Y = SPACE_normalOutput(Y_mask, X_mask, ROI_mask, pixel_size);
        else
            Results_Y = SPACE_nullOutput(Y_mask, ROI_mask);
        end
    end

end

function [Results] = SPACE_normalOutput(X_mask, Y_mask, ROI_mask, pixel_size)
% Characterize the spatial relationship of X relative to Y (X-->Y)

    Results = table;
    dt = bwdist(Y_mask) * pixel_size; % distance transformation of Y's mask
    Observed_NN_Distances = dt(X_mask & ROI_mask); % calculate observed nearest neighbor (NN) distances
    Random_NN_Distances = dt(ROI_mask); % calculate random nearest neighbor (NN) distances
    Results.Observed_Event_Count = numel(Observed_NN_Distances); % number of observed X events
    Results.Random_Event_Count = numel(Random_NN_Distances); % number of random X events
    [Results.Observed_x{1}, Results.Observed_PDF_y{1}] = epdf(Observed_NN_Distances); % compute observed PDF
    [Results.Random_x{1}, Results.Random_PDF_y{1}] = epdf(Random_NN_Distances); % compute random PDF
    Results.Observed_CDF_y{1} = pdf2cdf(Results.Observed_PDF_y{1}); % compute observed CDF
    Results.Random_CDF_y{1} = pdf2cdf(Results.Random_PDF_y{1}); % compute random CDF
    [Results.Delta_CDF_x{1}, Results.Delta_CDF_y{1}, Results.Spatial_Association_Index, Results.Spatial_Association_pValue] = CSRtest_single(Results.Observed_x{1}, Results.Observed_PDF_y{1}, Results.Observed_Event_Count, Results.Random_x{1}, Results.Random_PDF_y{1}, Results.Random_Event_Count); % compare observed and random CDFs

end

function [Results] = SPACE_nullOutput(X_mask, ROI_mask)
% At least one mask has no signal. Return null-result placeholder values
% for all table fields instead. Blank masks can be a common scenario, so
% this function provides an output that prevents errors and facilitates
% table concatenation. 

    Results = table;
    Results.Observed_Event_Count = sum(X_mask & ROI_mask,'all'); 
    Results.Random_Event_Count = sum(ROI_mask,'all');
    Results.Observed_x{1} = nan;
    Results.Observed_PDF_y{1} = nan;
    Results.Random_x{1} = nan;
    Results.Random_PDF_y{1} = nan;
    Results.Observed_CDF_y{1} = nan;
    Results.Random_CDF_y{1} = nan;
    Results.Delta_CDF_x{1} = nan;
    Results.Delta_CDF_y{1} = nan;
    Results.Spatial_Association_Index = nan;
    Results.Spatial_Association_pValue = nan;

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
    pdf_x = double(pdf_x);

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

function [delta_cdf_x, delta_cdf_y, spatial_association_index, p_value, verdict] = CSRtest_single(obs_pdf_x, obs_pdf_y, obs_n, ran_pdf_x, ran_pdf_y, ran_n, alpha)
% Test for complete spatial randomness (CSR) between 2 point patterns
% captured within a single image by comparing the observed
% (obs_pdf_x,obs_pdf_y) and random (ran_pdf_x,ran_pdf_y) distributions
% using their underlying PDF values. The p-value is calculated using their
% event sample sizes defined by obs_n and ran_n. Their delta function
% (delta_cdf_x,delta_cdf_y) is calulated as the random function subtraction
% from the observed function after they both are placed over a global
% x-coordinate scheme.
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/11/2023

    if nargin < 7
        alpha = 0.05; % default value
    end
    
    % Convert to column vector to ensure consistent output format
    obs_pdf_x = obs_pdf_x(:);
    obs_pdf_y = obs_pdf_y(:);
    ran_pdf_x = ran_pdf_x(:);
    ran_pdf_y = ran_pdf_y(:);
    
    % Regenerate input distributions over a global x-coordinate scheme
    delta_cdf_x = unique([obs_pdf_x; ran_pdf_x]); % find common bins
    [~, obs_idx] = ismember(obs_pdf_x, delta_cdf_x); % find where obs_pdf_y fits into global x-coordinates
    [~, ran_idx] = ismember(ran_pdf_x, delta_cdf_x); % find where ran_pdf_y fits into global x-coordinates
    obs_pdf_y_global = zeros(size(delta_cdf_x)); % initialize new y-values as zeros
    ran_pdf_y_global = obs_pdf_y_global; % initialize new y-values as zeros
    obs_pdf_y_global(obs_idx) = obs_pdf_y; % insert original y-values onto global x-coordinate scheme
    ran_pdf_y_global(ran_idx) = ran_pdf_y; % insert original y-values onto global x-coordinate scheme
    
    % Compute eCDFs from ePDFs
    obs_cdf_y = pdf2cdf(obs_pdf_y_global);
    ran_cdf_y = pdf2cdf(ran_pdf_y_global);
    
    % Compute the test statistic
    delta_cdf_y = obs_cdf_y - ran_cdf_y;
    [ks_p, ks_idx] = max(abs(delta_cdf_y));
    spatial_association_index = delta_cdf_y(ks_idx);
    
    % Compute the asymptotic P-value approximation and accept or reject the
    % null hypothesis on the basis of the P-value (based on code found in
    % built-in kstest2)
    n = obs_n * ran_n / (obs_n + ran_n);
    lambda = max((sqrt(n) + 0.12 + 0.11 / sqrt(n)) * ks_p, 0);
    p_value = exp(-2 * lambda * lambda);
    verdict = (alpha >= p_value);

end

function cdf_y = pdf2cdf(pdf_y)
% Compute an empirical cumulative distribution function (eCDF) based on
% an input empirical probability density function (ePDF).
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/11/2023

    cdf_y = cumsum(pdf_y);
    cdf_y = cdf_y ./ cdf_y(end);

end
