function [Results_Batch] = SPACE_batch(SPACE_Results)
% Perform global statistical analysis of a set of SPACE results. This will
% reproduce all functions over a globalized x-coordinate system such that
% each function will share the same x-coordinates without loss or
% modulation of the original y-coordinates, the function set's median and
% quantile envelope values, and the statistical significance of the set's
% spatial assocation.
% 
% -------------------------------------------------------------------------
%
% INPUTS:
%
% 1. SPACE_Results - SPACE results from individual analyses, formatted as a
%                    table with n-rows and 13-columns, where n is the
%                    number of individual SPACE analyses. 
%
% -------------------------------------------------------------------------
%
% OUTPUTS:
%
% 1. Results_Batch - Results summarizing the global statistical behavior of
%                    all individual SPACE analyses. Formatted as a table
%                    with 1-row and 21-columns. If all inidivual results
%                    are null, then batch analysis results will be
%                    undefined and most table fields will be set to 'NaN'.
%
% -------------------------------------------------------------------------
%
% AUTHORSHIP:
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/11/2023
%
% -------------------------------------------------------------------------


    % Input Validation  **************

    % Check for correct input data type
    if ~istable(SPACE_Results)
        error("Input is an invalid data type.");
    end

    % Check table dimensions
    nRows = height(SPACE_Results);
    nCols = width(SPACE_Results);
    if nRows < 1
        error("Input table is empty.");
    end
    if nCols < 13
        error("Input table has missing columns.");
    end
    if nCols > 13
        error("Input table has too many columns.");
    end

    % Check table fields
    if ~fieldsAreValid(SPACE_Results)
        error("Input table fields are invalid.");
    end


    % Perform Batch Analysis  **************

    % check if batch analysis is possible 
    batchAnalysisIsPossible = ~all(isnan(cell2mat(SPACE_Results.Observed_x)));

    Results_Batch = table;
    Results_Batch.Sample_Size = nRows;

    if batchAnalysisIsPossible    
        % Calculate global distributions
        Results_Batch.Global_x{1} = globalize_x(SPACE_Results.Observed_x, SPACE_Results.Random_x); % global x-coordinates
        [Results_Batch.Observed_CDF_x{1}, Results_Batch.Observed_CDF_y_Matrix{1}] = cdf_global(SPACE_Results.Observed_x, SPACE_Results.Observed_PDF_y, Results_Batch.Global_x{1}); % global observed CDF
        [Results_Batch.Random_CDF_x{1}, Results_Batch.Random_CDF_y_Matrix{1}] = cdf_global(SPACE_Results.Random_x, SPACE_Results.Random_PDF_y, Results_Batch.Global_x{1}); % global random CDF
        [Results_Batch.Delta_CDF_x{1}, Results_Batch.Delta_CDF_y_Matrix{1}] = delta_cdf_global(Results_Batch.Observed_CDF_y_Matrix{1}, Results_Batch.Random_CDF_y_Matrix{1}, Results_Batch.Global_x{1}); % global delta CDF
    
        % Calculate median functions and their quantile envelopes
        [Results_Batch.Observed_CDF_y_Median{1}, Results_Batch.Observed_CDF_y_Lower{1}, Results_Batch.Observed_CDF_y_Upper{1}] = median_fun(Results_Batch.Observed_CDF_y_Matrix{1}, Results_Batch.Observed_CDF_x{1}, SPACE_Results.Observed_Event_Count); % dispersed observed CDF
        [Results_Batch.Random_CDF_y_Median{1}, Results_Batch.Random_CDF_y_Lower{1}, Results_Batch.Random_CDF_y_Upper{1}] = median_fun(Results_Batch.Random_CDF_y_Matrix{1}, Results_Batch.Random_CDF_x{1}, SPACE_Results.Random_Event_Count); % dispersed random CDF
        [Results_Batch.Delta_CDF_y_Median{1}, Results_Batch.Delta_CDF_y_Lower{1}, Results_Batch.Delta_CDF_y_Upper{1}] = median_fun(Results_Batch.Delta_CDF_y_Matrix{1}, Results_Batch.Delta_CDF_x{1}, SPACE_Results.Observed_Event_Count); % dispersed delta CDF
    
        % Test point patterns for global CSR 
        Results_Batch.Sample_Event_Count{1} = SPACE_Results.Observed_Event_Count;
        Results_Batch.Sample_Spatial_Association_Index{1} = SPACE_Results.Spatial_Association_Index;
        [Results_Batch.Global_Spatial_Association_Index, Results_Batch.Global_CSR_Verdict] = CSRtest_batch(Results_Batch.Delta_CDF_y_Median{1}, Results_Batch.Delta_CDF_y_Lower{1}, Results_Batch.Delta_CDF_y_Upper{1}); % global CSR test
   
    else % results will be undefined, use placeholder values
        Results_Batch.Global_x{1} = nan;
        Results_Batch.Observed_CDF_x{1} = nan;
        Results_Batch.Observed_CDF_y_Matrix{1} = nan;
        Results_Batch.Random_CDF_x{1} = nan;
        Results_Batch.Random_CDF_y_Matrix{1} = nan;
        Results_Batch.Delta_CDF_x{1} = nan;
        Results_Batch.Delta_CDF_y_Matrix{1} = nan;
        Results_Batch.Observed_CDF_y_Median{1} = nan;
        Results_Batch.Observed_CDF_y_Lower{1} = nan;
        Results_Batch.Observed_CDF_y_Upper{1} = nan; 
        Results_Batch.Random_CDF_y_Median{1} = nan;
        Results_Batch.Random_CDF_y_Lower{1} = nan;
        Results_Batch.Random_CDF_y_Upper{1} = nan;
        Results_Batch.Delta_CDF_y_Median{1} = nan;
        Results_Batch.Delta_CDF_y_Lower{1} = nan;
        Results_Batch.Delta_CDF_y_Upper{1} = nan;
        Results_Batch.Sample_Event_Count{1} = nan;
        Results_Batch.Sample_Spatial_Association_Index{1} = nan;
        Results_Batch.Global_Spatial_Association_Index = nan;
        Results_Batch.Global_CSR_Verdict = nan;
    end

end

function [verdict] = fieldsAreValid(SPACE_Results)

    expectedFields = ["Observed_Event_Count", "Random_Event_Count",...
        "Observed_x", "Observed_PDF_y", "Random_x", "Random_PDF_y",...
        "Observed_CDF_y", "Random_CDF_y", "Delta_CDF_x", "Delta_CDF_y",...
        "Spatial_Association_Index", "Spatial_Association_pValue",...
        "Spatial_Association_Verdict"];

    actualFields = string(fields(SPACE_Results))';
    actualFields = actualFields(1:end-3); % ignore meta data fields 
    verdict = all(actualFields == expectedFields);

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
    x_global(isnan(x_global)) = []; % remove NaNs generated by masks with no signal, pdf data not useable
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

    % remove NaNs generated by masks with no signal, pdf data not useable
    badPDF = cellfun(@(pdf_x) all(isnan(pdf_x)), pdf_x_list);
    pdf_x_list(badPDF) = [];
    pdf_y_list(badPDF) = [];

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
% Calculate the weighted median and quantiles for an empricial function
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

    % THE FOLLOWING INTERMEDIATE VARIABLES WHICH ACT AS FUNCTION LISTS ARE
    % FORMATTED WITH ROWS CORRESPONDING TO EACH FUNCTION AND COLUMNS
    % CORRESPONDING TO FUNCTION Y-VALUES. THIS IS THE OPPOSITE FORMATTING
    % AS THE FUNCTION'S INPUT VARIABLES.

    % Sort function y-values for each x-coordinate
    [fun_y_sort_vals, fun_y_sort_idx] = sort(fun_y_mat', 1); % transpose for vectorized weighting later on
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
    weights_sum = cumsum(n_mat, 1); % cumulative sum of weights down columns
    weights_sum = weights_sum ./ sum(fun_n); % normalize to [0,1]

    % function median is first ordered function y-value whose
    % weight causes the weight-CDF to exceed or match 0.5
    median_thresh = 0.5;
    [~, median_row] = max(weights_sum >= median_thresh, [], 1); % find row subscripts of median values
    median_idx = median_row + rowNum2LinearIdx; % convert row numbers to linear indices
    fun_y_median = fun_y_sort_vals(median_idx)'; 

    % Give user option to set upper/lower as max/min
    if isnan(alpha) 
        fun_y_lower = min(fun_y_mat, [], 2);
        fun_y_upper = max(fun_y_mat, [], 2);
        return;
    end

    % Compute weighted quantiles, split alpha between lower and upper tails
    lower_thresh = alpha / 2; % compute quantile lower thresholds
    upper_thresh = 1 - lower_thresh; % compute quantile upper thresholds
    % Compute lower quantile
    [~, lower_row] = max(weights_sum >= lower_thresh, [], 1); % find row subscripts of lower quantile values
    lower_idx = lower_row + rowNum2LinearIdx; % convert row numbers to linear indices
    fun_y_lower = fun_y_sort_vals(lower_idx)'; % use linear indices to collect lower quantile values
    % Compute upper quantile
    [~, upper_row] = max(weights_sum >= upper_thresh, [], 1); % find row subscripts of upper quantile values
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
