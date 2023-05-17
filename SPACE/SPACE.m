function [varargout] = SPACE(X_mask, Y_mask, ROI_mask)
% Main function for performing 'Spatial Pattern Analysis using Closest
% Events' (SPACE) to characterize the spatial association between imaged
% patterns X and Y whose events are indicated by the 'true' pixels in their
% corresponding masks X_mask and Y_mask. All masks should be input as
% logical matrices. To analyze multiple images, input cell arrays of
% pattern masks. For optional subregion analysis, input a third mask which
% specifies a region-of-interest (ROI) within the image(s) indicated by the
% 'true' pixels in ROI_mask. If performing single image analysis, the
% results will be output as a single-row table with columns corresponding
% to the various single-image SPACE analysis results. If performing batch
% image analysis, the batch results will be the first output in the form of
% a single-row table with columns correspond to the varius SPACE batch
% analysis result, and the second output will be the single-image analysis
% results for each image in the batch in the form of a mutli-row table with
% rows correspind to each image and columns corresponding to the various
% single-image SPACE analysis results.
%
% See readme documentation for a detailed description of the output table
% fields.
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/11/2023

    % Check for appropriate number of inputs
    if nargin < 2
        error("Not enough input arguments.");
    elseif nargin > 3
        error("Too many input arguments.");
    end

    % Check if user is performing single-image analysis
    if islogical(X_mask) && islogical(Y_mask)
        if nargin == 2
            Single_Results = SPACE_single(X_mask, Y_mask, ROI_mask);
        elseif nargin == 3
            Single_Results = SPACE_single(X_mask, Y_mask, ROI_mask);
        end
        varargout = {Single_Results};
    end

    % Check if user is performing batch-image analysis
    if iscell(X_mask) && iscell(Y_mask)
        if nargin == 2
            [Batch_Results, Single_Results] = SPACE_batch(X_mask, Y_mask);
        elseif nargin == 3
            [Batch_Results, Single_Results] = SPACE_batch(X_mask, Y_mask, ROI_mask);
        end
        varargout = {Batch_Results, Single_Results};
    end

end

%% Single Image Analysis Functions

function [Results] = SPACE_single(X_mask, Y_mask, ROI_mask)
% Perform 'Spatial Pattern Analysis using Closest Events' (SPACE) of
% spatial point patterns X and Y whose events are identified by positive
% pixels in X_mask and Y_mask. Perform sub-region analysis by specifying a
% region of interest (ROI) identified by the positive pixels in ROI_mask. 
%
% See readme documentation for a detailed description of the output table
% fields.
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/11/2023

    % Input validation **************

    X_size = size(X_mask);
    Y_size = size(Y_mask);

    % Define default ROI
    if nargin < 3 % no ROI provided, create default
        ROI_mask = true(X_size);
    elseif nargin == 3
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

    % Compute distribution functions  **************

    Results = table;

    % Characterize the spatial relationship of X relative to Y (X-->Y)
    Y_dt = bwdist(Y_mask); % distance transformation of Y's mask
    XY_Observed_NN_Distances = Y_dt(X_mask & ROI_mask); % calculate observed nearest neighbor (NN) distances
    XY_Random_NN_Distances = Y_dt(ROI_mask); % calculate random nearest neighbor (NN) distances
    Results.XY_Observed_Event_Count = numel(XY_Observed_NN_Distances); % number of observed X events
    Results.XY_Random_Event_Count = numel(XY_Random_NN_Distances); % number of random X events
    [Results.XY_Observed_x{1}, Results.XY_Observed_PDF_y{1}] = epdf(XY_Observed_NN_Distances); % compute observed PDF
    [Results.XY_Random_x{1}, Results.XY_Random_PDF_y{1}] = epdf(XY_Random_NN_Distances); % compute random PDF
    Results.XY_Observed_CDF_y{1} = pdf2cdf(Results.XY_Observed_PDF_y{1}); % compute observed CDF
    Results.XY_Random_CDF_y{1} = pdf2cdf(Results.XY_Random_PDF_y{1}); % compute random CDF
    [Results.XY_Delta_CDF_x{1}, Results.XY_Delta_CDF_y{1}, Results.XY_Spatial_Association_Index, Results.XY_Spatial_Assocation_pValue] = CSRtest_single(Results.XY_Observed_x{1}, Results.XY_Observed_PDF_y{1}, Results.XY_Observed_Event_Count, Results.XY_Random_x{1}, Results.XY_Random_PDF_y{1}, Results.XY_Random_Event_Count); % compare observed and random CDFs

    % Characterize the spatial relationship of Y relative to X (Y-->X)
    X_dt = bwdist(X_mask); % distance transformation of X's mask
    YX_Observed_NN_Distances = X_dt(Y_mask & ROI_mask); % calculate observed nearest neighbor (NN) distances
    YX_Random_NN_Distances = X_dt(ROI_mask); % calculate random nearest neighbor (NN) distances
    Results.YX_Observed_Event_Count = numel(YX_Observed_NN_Distances); % number of observed Y events
    Results.YX_Random_Event_Count = numel(YX_Random_NN_Distances); % number of random y events
    [Results.YX_Observed_x{1}, Results.YX_Observed_PDF_y{1}] = epdf(YX_Observed_NN_Distances); % compute observed PDF
    [Results.YX_Random_x{1}, Results.YX_Random_PDF_y{1}] = epdf(YX_Random_NN_Distances); % compute random PDF
    Results.YX_Observed_CDF_y{1} = pdf2cdf(Results.YX_Observed_PDF_y{1}); % compute observed CDF
    Results.YX_Random_CDF_y{1} = pdf2cdf(Results.YX_Random_PDF_y{1}); % compute random CDF
    [Results.YX_Delta_CDF_x{1}, Results.YX_Delta_CDF_y{1}, Results.YX_Spatial_Association_Index, Results.YX_Spatial_Assocation_pValue] = CSRtest_single(Results.YX_Observed_x{1}, Results.YX_Observed_PDF_y{1}, Results.YX_Observed_Event_Count, Results.YX_Random_x{1}, Results.YX_Random_PDF_y{1}, Results.YX_Random_Event_Count); % compare observed and random CDFs

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

function [delta_cdf_x, delta_cdf_y, spatial_assocation_index, p_value, verdict] = CSRtest_single(obs_pdf_x, obs_pdf_y, obs_n, ran_pdf_x, ran_pdf_y, ran_n, alpha)
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
    obs_pdf_x =  obs_pdf_x(:);
    obs_pdf_y =  obs_pdf_y(:);
    ran_pdf_x  =  ran_pdf_x(:);
    ran_pdf_y  =  ran_pdf_y(:);

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
    spatial_assocation_index = delta_cdf_y(ks_idx);
    
    % Compute the asymptotic P-value approximation and accept or reject the
    % null hypothesis on the basis of the P-value (based on code found in
    % built-in kstest2)
    n = obs_n * ran_n / (obs_n + ran_n);
    lambda = max((sqrt(n) + 0.12 + 0.11 / sqrt(n)) * ks_p, 0);
    p_value  =  exp(-2 * lambda * lambda);
    verdict = (alpha >= p_value);

end

%% Batch Image Analysis Functions

function [Batch_Results, Single_Results] = SPACE_batch(X_mask_list, Y_mask_list, ROI_mask_list)
% Perform 'Spatial Pattern Analysis using Closest Events' (SPACE) on
% multiple samples of spatial point patterns X and Y whose events are
% identified by positive pixels in the cell array of image masks,
% X_mask_list and Y_mask_list. Perform sub-region analysis by specifying a
% region of interest (ROI) identified by the positive pixels in each mask
% contained in the cell array ROI_mask_list. Batch results are output in
% the SPACE_Group table and individual image results are output in the
% SPACE_Sample table. 
%
% See readme documentation for a detailed description of the output table
% fields.
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/11/2023

    % Input validation  **************

    % Convert everything to column vectors for consistent formatting
    X_mask_list = X_mask_list(:);
    Y_mask_list = Y_mask_list(:);
    ROI_mask_list = ROI_mask_list(:);

    X_mask_count = numel(X_mask_list);
    Y_mask_count = numel(Y_mask_list);
    ROI_mask_count = numel(ROI_mask_list);

    % Make sure there are the same number of masks in each list
    if numel(unique([X_mask_count, Y_mask_count, ROI_mask_count])) > 1
        error("Mask cell arrays are of different lengths.");
    end
    image_count = X_mask_count;

    % Make sure each mask is a cell array
    if ~iscell(X_mask_list) || ~ iscell(Y_mask_list) || ~ iscell(ROI_mask_list)
        error("Mask lists must be cell arrays");
    end

    % Perform single sample SPACE   **************
    
    Single_Results = table;
    wb = waitbar(0, sprintf("Analyzing image: %i/%i", 0, image_count));
    for i_image = 1:image_count
        Single_Results(i_image,:) = SPACE_single(X_mask_list{i_image}, Y_mask_list{i_image}, ROI_mask_list{i_image});
        waitbar(i_image/image_count, wb, sprintf("Analyzing image: %i/%i", i_image, image_count));
    end
    close(wb);

    Batch_Results = table;
    Batch_Results.Sample_Size = image_count;

    % X-->Y Batch Analysis  **************

    % Calculate global distributions
    Batch_Results.XY_Global_x{1} = globalize_x(Single_Results.XY_Observed_x, Single_Results.XY_Random_x); % global x-coordinates
    [Batch_Results.XY_Observed_CDF_x{1}, Batch_Results.XY_Observed_CDF_y_Matrix{1}] = cdf_global(Single_Results.XY_Observed_x, Single_Results.XY_Observed_PDF_y, Batch_Results.XY_Global_x{1}); % global observed CDF
    [Batch_Results.XY_Random_CDF_x{1}, Batch_Results.XY_Random_CDF_y_Matrix{1}] = cdf_global(Single_Results.XY_Random_x, Single_Results.XY_Random_PDF_y, Batch_Results.XY_Global_x{1}); % global random CDF
    [Batch_Results.XY_Delta_CDF_x{1}, Batch_Results.XY_Delta_CDF_y_Matrix{1}] = delta_cdf_global(Batch_Results.XY_Observed_CDF_y_Matrix{1}, Batch_Results.XY_Random_CDF_y_Matrix{1}, Batch_Results.XY_Global_x{1}); % global delta CDF

    % Calculate median functions and their quantile envelopes
    [Batch_Results.XY_Observed_CDF_y_Median{1}, Batch_Results.XY_Observed_CDF_y_Lower{1}, Batch_Results.XY_Observed_CDF_y_Upper{1}] = median_fun(Batch_Results.XY_Observed_CDF_y_Matrix{1}, Batch_Results.XY_Observed_CDF_x{1}, Single_Results.XY_Observed_Event_Count); % dispersed observed CDF
    [Batch_Results.XY_Random_CDF_y_Median{1}, Batch_Results.XY_Random_CDF_y_Lower{1}, Batch_Results.XY_Random_CDF_y_Upper{1}] = median_fun(Batch_Results.XY_Random_CDF_y_Matrix{1}, Batch_Results.XY_Random_CDF_x{1}, Single_Results.XY_Random_Event_Count); % dispersed random CDF
    [Batch_Results.XY_Delta_CDF_y_Median{1}, Batch_Results.XY_Delta_CDF_y_Lower{1}, Batch_Results.XY_Delta_CDF_y_Upper{1}] = median_fun(Batch_Results.XY_Delta_CDF_y_Matrix{1}, Batch_Results.XY_Delta_CDF_x{1}, Single_Results.XY_Observed_Event_Count); % dispersed delta CDF

    % Test point patterns for global CSR 
    Batch_Results.XY_Sample_Event_Count{1} = Single_Results.XY_Observed_Event_Count;
    Batch_Results.XY_Sample_Spatial_Association_Index{1} = Single_Results.XY_Spatial_Association_Index;
    [Batch_Results.XY_Global_Spatial_Association_Index, Batch_Results.XY_Global_CSR_Verdict] = CSRtest_batch(Batch_Results.XY_Delta_CDF_y_Median{1}, Batch_Results.XY_Delta_CDF_y_Lower{1}, Batch_Results.XY_Delta_CDF_y_Upper{1}); % global CSR test

    % Y-->X Batch Analysis  **************

    % Calculate global distributions
    Batch_Results.YX_Global_x{1} = globalize_x(Single_Results.YX_Observed_x, Single_Results.YX_Random_x); % global x-coordinates
    [Batch_Results.YX_Observed_CDF_x{1}, Batch_Results.YX_Observed_CDF_y_Matrix{1}] = cdf_global(Single_Results.YX_Observed_x, Single_Results.YX_Observed_PDF_y, Batch_Results.YX_Global_x{1}); % global observed CDF
    [Batch_Results.YX_Random_CDF_x{1}, Batch_Results.YX_Random_CDF_y_Matrix{1}] = cdf_global(Single_Results.YX_Random_x, Single_Results.YX_Random_PDF_y, Batch_Results.YX_Global_x{1}); % global random CDF
    [Batch_Results.YX_Delta_CDF_x{1}, Batch_Results.YX_Delta_CDF_y_Matrix{1}] = delta_cdf_global(Batch_Results.YX_Observed_CDF_y_Matrix{1}, Batch_Results.YX_Random_CDF_y_Matrix{1}, Batch_Results.YX_Global_x{1}); % global delta CDF

    % Calculate median functions and their quantile envelopes
    [Batch_Results.YX_Observed_CDF_y_Median{1}, Batch_Results.YX_Observed_CDF_y_Lower{1}, Batch_Results.YX_Observed_CDF_y_Upper{1}] = median_fun(Batch_Results.YX_Observed_CDF_y_Matrix{1}, Batch_Results.YX_Observed_CDF_x{1}, Single_Results.YX_Observed_Event_Count); % dispersed observed CDF
    [Batch_Results.YX_Random_CDF_y_Median{1}, Batch_Results.YX_Random_CDF_y_Lower{1}, Batch_Results.YX_Random_CDF_y_Upper{1}] = median_fun(Batch_Results.YX_Random_CDF_y_Matrix{1}, Batch_Results.YX_Random_CDF_x{1}, Single_Results.YX_Random_Event_Count); % dispersed random CDF
    [Batch_Results.YX_Delta_CDF_y_Median{1}, Batch_Results.YX_Delta_CDF_y_Lower{1}, Batch_Results.YX_Delta_CDF_y_Upper{1}] = median_fun(Batch_Results.YX_Delta_CDF_y_Matrix{1}, Batch_Results.YX_Delta_CDF_x{1}, Single_Results.YX_Observed_Event_Count); % dispersed delta CDF

    % Test point patterns for global CSR 
    Batch_Results.YX_Sample_Event_Count{1} = Single_Results.YX_Observed_Event_Count;
    Batch_Results.YX_Sample_Spatial_Association_Index{1} = Single_Results.YX_Spatial_Association_Index;
    [Batch_Results.YX_Global_Spatial_Association_Index, Batch_Results.YX_CSR_Verdict] = CSRtest_batch(Batch_Results.YX_Delta_CDF_y_Median{1}, Batch_Results.YX_Delta_CDF_y_Lower{1}, Batch_Results.YX_Delta_CDF_y_Upper{1}); % global CSR test

end

function x_global = globalize_x(varargin)
% Find unique x-coordinates and sort in ascending order. Input
% comma-separated list of x-coordinates for any number of epirical
% functions.
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/11/2023

    varargin = varargin(:); % convert to column vector for consistent formatting
    x_global = cell2mat(cellfun(@(x) cell2mat(x), varargin, 'uniformoutput', false)); % unwind the arrays
    x_global = unique(x_global); % return ordered unique x-coordinates

end

function [cdf_x, cdf_y_mat] = cdf_global(pdf_x_list, pdf_y_list, x_global)
% Calculate the global cumulate distribution functions (CDFs) for a list of
% individual functions defined over a global x-coordinate scheme (x_global)
% such that each function is defined over the same x-coordinates.
% Individual dsitribution function coordinates (pdf_x_list and pdf_y_list)
% should be stored in cell array of column vectors. The output global CDFs
% are stored in a matrix with individual functions stored in each column
% and their y-coordinates stored across rows; the global x-coordinates for
% all functions are stored the array cdf_x. For example, the global x- and
% y-coordinates for function 'i' can be found at cdf_x(:) and
% cdf_y_mat(:,i).
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/11/2023

    % Convert to column vector to ensure consistent output format
    pdf_x_list = pdf_x_list(:);
    pdf_y_list = pdf_y_list(:);
    x_global = x_global(:);
    pdf_x_list = cellfun(@(x) x(:), pdf_x_list, 'uniformoutput', false);
    pdf_y_list = cellfun(@(y) y(:), pdf_y_list, 'uniformoutput', false);

    % count number of samples
    n_samples = numel(pdf_x_list);

    % Get common x-coordinates
    cdf_x = unique(cell2mat(pdf_x_list)); % find all unique x-coordinates between all samples
    if nargin > 2
        if ~isempty(x_global)
            % if a global x-coordinate scheme is provided, use only the relevant values
            x_max = cdf_x(end);
            cdf_x = x_global(x_global <= x_max);
        end
    end
    n_x = numel(cdf_x); % count the total number of x-coordinates

    % Put PDFs on global x_coordinates and convert to CDF
    pdf_y_blank = zeros(n_x, 1); % preallocate global y-range once and reuse each loop
    cdf_y_mat = nan(n_x, n_samples); % preallocate matrix to store all CDF y-coordinates
    for i_sample = 1:n_samples
        pdf_x = pdf_x_list{i_sample}(:); % make sure this is a column vector
        pdf_y = pdf_y_list{i_sample}(:); % make sure this is a column vector
        [~,x_idx] = ismember(pdf_x, cdf_x); % find the x-values for this distribution's y-coordinates
        pdf_y_common = pdf_y_blank; % intantiate PDF y-coordinates as all zeros
        pdf_y_common(x_idx) = pdf_y; % insert this sample's PDF y-coordinates
        cdf_y_mat(:, i_sample) = pdf2cdf(pdf_y_common); % compute this sample's CDF and add it to the global matrix
    end

end

function [delta_cdf_x, delta_cdf_y_mat] = delta_cdf_global(obs_cdf_mat, ran_cdf_mat, x_global)
% Calculate the global delta functions for paired samples of observed and
% random distribution functions. Each set of functions should be input as
% matrices (obs_cdf_mat and ran_cdf_mat) of functions defined over the same
% global x-coordinate scheme (x_global). The rows of these matrices should
% correspond to function y-coordinates while columns should correspond to
% individual function. For example, the x- and y-coordinates of observed
% function 'i' can be found at x_global(:) and obs_cdf_mat(:,i). Individual
% delta functions are defined over a global x-coordinate scheme
% (delta_cdf_x) and y-coordinates for individual functions are also stored
% in a matrix with rows for y-coordinates and columns for individual
% functions. For example, the delta function x- and y-coordinates for
% sample 'i' can be found at delta_cdf_x(:) and delta_cdf_y_mat(:,i). 
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/11/2023

    % Convert to column vector to ensure consistent output format
    x_global = x_global(:);

    % Input Valdiation
    [n_x_obs,n_samples_obs] = size(obs_cdf_mat);
    [n_x_ran,n_samples_ran] = size(ran_cdf_mat);
    if n_samples_obs ~= n_samples_ran
        error("Observed and random sample sizes do not match.");
    end

    % get everything on the same x-range
    n_samples = n_samples_obs;
    y_padding = ones(abs(n_x_ran - n_x_obs), n_samples); % reuse padding
    if n_x_obs > n_x_ran
        ran_cdf_mat = [ran_cdf_mat; y_padding];
        delta_cdf_x = x_global(1:n_x_obs);
    elseif n_x_obs < n_x_ran
        obs_cdf_mat = [obs_cdf_mat; y_padding];
        delta_cdf_x = x_global(1:n_x_ran);
    else
        delta_cdf_x = x_global(1:n_x_obs);
    end

    % Compare observed and random distributions
    delta_cdf_y_mat = obs_cdf_mat - ran_cdf_mat;

end

function [fun_y_median, fun_y_lower, fun_y_upper] = median_fun(fun_y_mat, fun_x, fun_n, alpha)
% Calculated the weighted median and quantiles for an empricial function
% based on individual functions defined by fun_y_mat, where columns
% correspond to each function and rows correspond to function y-values. All
% individual functions captured within fun_y_mat are assumed to be defined
% over the same x-coordinates, captured by fun_x. Individual function
% weights are calculated as the function's sample size, captured by fun_n,
% divided by the total sample size across all individual functions (ie
% the weight of individual function 'i' is weight(i) = fun_n(i) /
% sum(fun_n). If no alpha value is specified, a default value of 0.05 will
% be used.
%
% Weighted median and quantiles are calculated according to the Wikipedia
% page - https://en.wikipedia.org/wiki/Weighted_median. 
% sample 'i' can be found at delta_cdf_x(:) and delta_cdf_y_mat(:,i). 
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/11/2023

    % Input validation
    if nargin < 4
        alpha = 0.1;  % default alpha value is 10%, split to 5% for each tail
    else
        if (alpha > 1 || alpha < 0) || ~isscalar(alpha)
            error("Alpha must be a scalar number between 0 and 1, inclusive.");
        end
    end

    % get everything on the same x-range
    x_count = numel(fun_x);
    fun_count = size(fun_y_mat, 2);

    % give functions equal weight if only one weight is supplied
    n_count = numel(fun_n);
    if (fun_count ~= n_count) && (n_count == 1)
        fun_n = repelem(fun_n ,fun_count, 1);
    end

    % Sort function y-values for each x-coordinate
    [fun_y_sort_vals, fun_y_sort_idx] = sort(fun_y_mat'); % transpose for vectorized weighting later on
    rowNum2LinearIdx = (0:(x_count-1)) * fun_count; % indices returned as row numbers, needs conversion to linear index by adding idx*colNumber to row numbers
    fun_y_sort_idx = fun_y_sort_idx + repmat(rowNum2LinearIdx, fun_count, 1); % convert row numbers to linear indices

    % Repeat individual sample sizes for each x-coordinate and sort
    % according to function y-value scheme from above
    fun_n = fun_n(:); % make sure this is a column vector
    n_mat = repmat(fun_n, 1, x_count); % row of sample sizes per x-coordinate
    n_mat = n_mat(fun_y_sort_idx); % order weights based on CDF y-value sorting

    % Calculte weight of each function as each sample's observation
    % count divided by the total observation count (sum of all observation
    % counts).
    weights_sum = cumsum(n_mat); % cumulative sum of weights down rows
    weights_sum = weights_sum ./ sum(fun_n); % normalize to [0,1]

    % function median is first ordered function y-value whose
    % weight causes the weight-CDF to exceed or match 0.5
    median_thresh = 0.5;
    [~, median_row] = max(weights_sum >= median_thresh); % find row numbers of median values
    median_idx = median_row + rowNum2LinearIdx; % convert row numbers to linear indices
    fun_y_median = fun_y_sort_vals(median_idx)'; 

    if isnan(alpha) % give user option to set upper/lower as max/min
        fun_y_lower = min(fun_y_mat, [], 2);
        fun_y_upper = max(fun_y_mat, [], 2);
        return;
    end

    % Compute weighted quantiles, split alpha between lower and upper tails
    % Compute quantile thresholds
    lower_thresh = alpha / 2;
    upper_thresh = 1 - lower_thresh;
    % Compute lower quantile
    [~, lower_row] = max(weights_sum >= lower_thresh); % find row numbers of lower quantile values
    lower_idx = lower_row + rowNum2LinearIdx; % convert row numbers to linear indices
    fun_y_lower = fun_y_sort_vals(lower_idx)'; % use linear indices to collect lower quantile values
    % Compute upper quantile
    [~, upper_row] = max(weights_sum >= upper_thresh); % find row numbers of upper quantile values
    upper_idx = upper_row + rowNum2LinearIdx; % convert row numbers to linear indices
    fun_y_upper = fun_y_sort_vals(upper_idx)'; % use linear indices to collect upper quantile values

    % Plot results for troubleshooting and validation 
    % figure;
    % hold on;
    % c = plot(fun_x,fun_y_mat,'r-');
    % p(1) = plot(fun_x, fun_y_median, 'k', 'linewidth', 4);
    % p(2) = plot(fun_x, fun_y_lower, 'k--', 'linewidth', 3);
    % p(3) = plot(fun_x, fun_y_upper, 'k:', 'linewidth', 3);
    % legend([c(1), p], ["Individual", "Median", "Lower", "Upper"]);
    % grid on;
    % hold off; 

end

function [spatial_assocation_index, verdict] = CSRtest_batch(delta_cdf_median, delta_cdf_lower, delta_cdf_upper)
% Test for complete spatial randomness (CSR) between 2 point patterns
% captured by multiple images by evaluating whether the qunatile envelope
% of the median delta function overlaps with 0. If the envelope psoitively
% exceeds 0, there is sufficient evidence to conclude the point patterns
% are aggregated. If the envelope negatively exceeds 0, then there is
% sufficient evidence to conclude that the point patterns are dispersed.
% And if the envelope overlaps 0 for the entire distance range, then there
% is not sufficient evidence to differentiate the point patterns'
% relationship from CSR. If there is deviation from zero, the global
% spatial association index is defined as the value of the median function
% at the absolute maximal deviation of the quantile envelope. Otherwise,
% global spatial association is defined as the absolute maximum of the
% median function.
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/11/2023

    lower_deviation = delta_cdf_lower > 0;
    upper_deviation = delta_cdf_upper < 0;
    
    if any(lower_deviation) || any(upper_deviation)
        verdict = true; % deviation from 0 detected

        % Calculate spatial association index
        [lower_deviation_maxVal, lower_deviation_maxIdx] = max(abs(delta_cdf_lower(lower_deviation)));
        [upper_deviation_maxVal, upper_deviation_maxIdx] = max(abs(delta_cdf_upper(upper_deviation)));
        if isempty(lower_deviation_maxVal)
            % no deviation here causes empty matrix and inequaily fails
            lower_deviation_maxVal = 0;
        end
        if isempty(upper_deviation_maxVal)
            % no deviation here causes empty matrix and inequaily fails
            upper_deviation_maxVal = 0;
        end
        if lower_deviation_maxVal > upper_deviation_maxVal
            spatial_assocation_index = delta_cdf_median(lower_deviation_maxIdx);
        else
            spatial_assocation_index = delta_cdf_median(upper_deviation_maxIdx);
        end
    else
        verdict = false; % no deviation from 0 
        spatial_assocation_index = max(delta_cdf_median, [], 'comparisonmethod', 'abs');
    end 

end

%% Shared Functions

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
