function [Results] = SPACE(mask_subject, mask_landmark, mask_ROI, pixelSize, alpha)
% Perform 'Spatial Pattern Analysis using Closest Events' (SPACE), a
% discrete, bivariate, nearest neighbor-based point patten analysis which
% characterizes the spatial association between two different patterns
% within a digital image. Analysis results describe the association of the
% first (subject) pattern's events relative to the second (landmark)
% pattern's events.
%
% 
% -------------------------------------------------------------------------
%
% INPUTS:
%
% 1. MASK_SUBJECT - Discrete events which compose an imaged pattern whose
%                   spatial distribution relative to the landmark pattern
%                   will be described by the analysis results. Formatted as
%                   a logical matrix (binary image mask) of any dimension
%                   with events represented by 'true' elements. The matrix
%                   size must match the sizes of MASK_LANDMARK and
%                   MASK_ROI.
%
% 2. MASK_LANDMARK - Discrete events which compose an imaged pattern that
%                    will serve as the landmark by which the subject
%                    pattern's spatial distribution will be derived.
%                    Formatted as a logical matrix (binary image mask) of
%                    any dimension with events represented by 'true'
%                    elements. The matrix size must match the sizes of
%                    MASK_SUBJECT and MASK_ROI.
%
% 3. MASK_ROI - OPTIONAL input which is a logical matrix (binary image
%               mask) of any dimension used to indicate a
%               region-of-interest (ROI) within the full study area where
%               the analysis will be exclusively performed. Event pixels
%               which fall outside of this region will contribute no
%               information to the analysis. If you do not want to specify
%               an ROI but do want to provide subsequent function inputs,
%               then provide this input as an empty matrix. If this input
%               is not provided or input as an empty matrix, then the
%               entire study area will be considered and all event pixels
%               will be used in the analysis. Here, the set of elements
%               with a value of 1 or 'true' represent the region(s) to
%               include in the analysis and the set of elements with a
%               value of 0 or 'false' represent the region(s) to omit from
%               the analysis. This matrix must have the same size as
%               MASK_SUBJECT and MASK_LANDMARK. This input is MANDATORY if
%               MASK_SUBJECT and MASK_LANDMARK were first processed using
%               the Isotropic Replacement method, otherwise this analyses
%               will produce 100% garbage results.
%
% 4. PIXELSIZE - OPTIONAL input which specifies the real-world size
%                (side-length) of pixels composing the digital image(s)
%                being analyzed. Formatted as a positive numeric scalar. If
%                you do not want to specify a pixel size but do want to
%                provide subsequent function inputs, then provide this
%                input as an empty matrix. If this input is not provided or
%                input as an empty matrix, then distance-derived outputs
%                will have units of 'pixel-count'. This input is used to
%                scale SPACE results so they have the same-real world units
%                as the input digital image(s). For the SPACE analysis to
%                produce results which are not 100% garbage, the real-world
%                shape of image pixels must be cuboidal, meaning every side
%                of the pixel has the same length, and this length is the
%                single value that needs to be assigned to this input
%                variable. If your pixels are in fact not cuboidal (e.g.
%                the x/y-resolution of your digital image is different than
%                its z-resolution), fear not! You just have to resample the
%                masks so that they are isotropic and have cuboidal pixels.
%                This can easily be done using the Isotropic Replacement
%                method found at the GitHub link below. From this method,
%                the 'isotropic_ROI mask' (2nd) output must then be used as
%                the MASK_ROI input to this function for the SPACE analysis
%                to be performed correctly.
%
% 5. ALPHA - OPTIONAL input which specifies the significance level for
%            testing complete spatial randomness of the subject pattern
%            relative to the landmark pattern. Formatted as a numeric
%            scalar with a value between 0 and 1, exlusive. A value of 0.05
%            is used by default if no value is provided.
%
%
% -------------------------------------------------------------------------
%
% OUTPUTS: 
%
% 1. RESULST - SPACE analysis results formatted as a 1-by-12 table which
%              describe the spatial distribution of the subject pattern's
%              events relative to the landmark pattern's events. See the
%              linked GitHub page below for info on what each table field
%              represents. If at least one of the patterns has no events,
%              then the analysis cannot be performed and most table fields
%              will automatically be assigned a value of 'NaN' to avoid
%              throwing errors and crashing batch analyses. NOTE: To fully
%              describe the spatial association between two different
%              patterns, this analyis must be performed twice, such that
%              each pattern has been treated as the subject and landmark
%              pattern.
%
%
% -------------------------------------------------------------------------
%
% USES CASES:
% 
% [Results] = SPACE(mask_subject, mask_landmark) 
% Perform SPACE using default ROI, pixel size, and alpha. The entire study
% space will be evaluated, pixels will be assigned a size (side length) of
% 1, and a significance level (alpha) of 0.05 will be used to evaluate
% complete spatial randomness.
%
% [Results] = SPACE(mask_subject, mask_landmark, ROI_mask) 
% Perform SPACE on a subset of the full study area as specified by
% MASK_ROI. Pixels are assigned a default size of 1.
%
% [Results] = SPACE(mask_subject, mask_landmark, [], pixelSize) 
% Perform SPACE on the full study area. Pixels are assigned a size
% according to PIXELSIZE and distance-based output results are scaled
% accordingly.
%
% [Results] = SPACE(mask_subject, mask_landmark, mask_ROI, pixelSize) 
% Perform SPACE on a user-defined subset of the full study area and pixel
% size.
%
% [Results] = SPACE(mask_subject, mask_landmark, ___, ___, alpha) 
% Perform SPACE with a user-defined significance level for determining
% complete spatial randomess. 
%
%
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
% Last Updated: 03/27/2025
%
% Copyright 2025, Andrew Michael Soltisz. All rights reserved.
% This source code is licensed under the BSD-3-Clause License found in the
% LICENSE.txt file in the root directory of this source tree.
%
%
% -------------------------------------------------------------------------


    % --------------- INPUT VALIDATION ---------------

    % check for appropriate number of inputs
    if nargin < 2
        error("Not enough input arguments.");
    end

    % validate MASK_SUBJECT and MASK_LANDMARK inputs
    if ~islogical(mask_subject)
        error("Input MASK_SUBJECT must be logical matrix.");
    end
    if ~islogical(mask_landmark)
        error("Input MASK_LANDMARK must be logical matrix.");
    end

    % validate image sizes
    imageSize_subject = size(mask_subject);
    imageSize_landmark = size(mask_landmark);
    if any(imageSize_subject ~= imageSize_landmark)
        error("Inputs MASK_SUBJECT and MASK_LANDMARK must have the same size.");
    end

    % validate input MASK_ROI
    if nargin < 3 % no ROI provided, use default
        mask_ROI = true(imageSize_subject);
    elseif nargin >= 3
        if isempty(mask_ROI) % ROI is empty, create default
            mask_ROI = true(imageSize_subject);
        else
            if ~islogical(mask_ROI)
                error("Input MASK_ROI must be a logical matrix.");
            end
            if ~all(imageSize_subject == size(mask_ROI))
                error("Input MASK_ROI size does not match image size of inputs MASK_SUBJECT or MASK_LANDMARK.");
            end
        end
    end

    % validate input PIXELSIZE
    if nargin < 4 % no pixel size provided, use default
        pixelSize = 1; % default value
    else
        if isempty(pixelSize) % pixel size is empty, use default
            pixelSize = 1; % default value
        else
            if ~(isscalar(pixelSize) && isnumeric(pixelSize))
                error("Input PIXELSIZE must be a numeric scalar.");
            end
            if pixelSize == 0
                error("Input PIXELSIZE must be non-zero.");
            end
        end
    end

    % validate input ALPHA
    if nargin < 5 % no alpha provided, use default
        alpha = 0.05; % default value
    else
        if ~(isscalar(alpha) && isnumeric(alpha))
            error("Input ALPHA must be a numeric scalar.");
        end
        if alpha <= 0 || alpha >= 1
            error("Input ALPHA value must be between 0 and 1 (non-inclusive).");
        end
    end


    % --------------- ANALYSIS ---------------

    % remove out-of-bounds-pixels
    mask_subject = mask_subject & mask_ROI;
    mask_landmark = mask_landmark & mask_ROI;
    noNullInputs = any(mask_subject,'all') && any(mask_landmark,'all'); % check for empty masks

    Results = table; % create output table data structure
    if noNullInputs
        % measure inter-event distances
        distanceTransform_landmark = bwdist(mask_landmark) .* pixelSize; % distance transformation of landmark mask
        observedDistances = distance_pixel2pixel(mask_subject, distanceTransform_landmark); % calculate observed nearest neighbor (NN) distances
        randomDistances = distance_pixel2pixel(mask_ROI, distanceTransform_landmark); % calculate random nearest neighbor (NN) distances
    
        % count events
        Results.Observed_Event_Count = numel(observedDistances); % count number of observed subject events (pixels)
        Results.Random_Event_Count = numel(randomDistances); % count number of random subject events (pixels)
    
        % compute empirical probability density functions (ePDFs)
        [Results.Observed_x{1}, Results.Observed_PDF_y{1}] = ePDF(observedDistances); % compute observed PDF coordintes
        [Results.Random_x{1}, Results.Random_PDF_y{1}] = ePDF(randomDistances); % compute random PDF coordintes
        
        % compute cumulative distribution functions (CDFs)
        Results.Observed_CDF_y{1} = PDF2CDF(Results.Observed_PDF_y{1}); % compute observed CDF coordintes
        Results.Random_CDF_y{1} = PDF2CDF(Results.Random_PDF_y{1}); % compute random CDF coordintes
    
        % test for complete spatial randomness (CSR) by comparoing observed and random CDFs
        [Results.Delta_CDF_x{1}, Results.Delta_CDF_y{1}, Results.Spatial_Association_Index, Results.Spatial_Association_pValue, Results.Spatial_Association_Verdict] = testCSR_pixel(Results.Observed_x{1}, Results.Observed_PDF_y{1}, Results.Observed_Event_Count, Results.Random_x{1}, Results.Random_PDF_y{1}, Results.Random_Event_Count, alpha); % compare observed and random CDFs
    else
        % At least one species has no events. Return null-result
        % placeholder values for all table fields. This functionality
        % provides an output to prevent errors from halting batch analyses
        % while facilitating concatenation of results table rows.
        Results.Observed_Event_Count = sum(mask_subject & mask_landmark, 'all'); 
        Results.Random_Event_Count = sum(mask_landmark, 'all');
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
        Results.Spatial_Association_Verdict = nan;
    end

end

%% Helper Functions

function [distances_subjectPixels] = distance_pixel2pixel(mask_subject, distanceTransform_landmark)
% Measure the distance from each subject pixel to its nearest neighbor (NN)
% landmark pixel.
%
% 
% -------------------------------------------------------------------------
%
% INPUTS:
%
% 1. MASK_SUBJECT - Logical matrix (binary image mask) of any dimension
%                   which defines the subject pattern with events
%                   represented by 'true' elements. This matrix must have
%                   the same size as DISTANCETRANSFORM_LANDMARK.
%
% 2. DISTANCETRANSFORM_LANDMARK - Numeric matrix that is the distance
%                                 transformation of the landmark pattern's
%                                 binary image mask. The value of each
%                                 element specifies the Euclidean distance
%                                 of that pixel to the nearest event pixel
%                                 in the landmark pattern's mask.
%
%
% -------------------------------------------------------------------------
%
% OUTPUTS: 
%
% 1. DISTANCES_SUBJECTPIXELS - Numeric column vector listing the distance
%                              of each 'true' pixel in MASK_SUBJECT to the
%                              nearest 'true' pixel in the landmark
%                              pattern's mask. This vector will have as
%                              many elements as there are 'true' pixels in
%                              MASK_SUBJECT.
%
%
% -------------------------------------------------------------------------

    distances_subjectPixels = distanceTransform_landmark(mask_subject); % calculate observed nearest neighbor (NN) distances
    distances_subjectPixels = distances_subjectPixels(:); % ensure column fector format

end

function [delta_CDF_x, delta_CDF_y, testStatistic, pValue, verdict] = testCSR_pixel(observed_PDF_x, observed_PDF_y, observed_nEvents, random_PDF_x, random_PDF_y, random_nEvents, alpha)
% Test for complete spatial randomness (CSR) between two point patterns
% captured within a single image using a 2-sided Kolmogorov-Smirnov test
% which evalutes equality of the cumulative distribution functions of
% observed and random distances. The pValue is calculated using event
% sample sizes. The function describing the difference between these CDFs,
% called the delta CDF , is calulated as the random CDF subtracted from the
% observed CDF after they both are placed over a global x-coordinate
% scheme.
%
% 
% -------------------------------------------------------------------------
%
% INPUTS:
%
% 1. OBSERVED_PDF_X - Numeric column vector listing the x-coordinates
%                     (representing distance) of the PDF of the observed
%                     nearest neighbor distances between each subject event
%                     and their nearest landmark event. This must have as
%                     many elements as the input OBSERVED_PDF_Y and The
%                     elemenets must monotonically increase.
%
% 2. OBSERVED_PDF_Y - Numeric column vector listing the x-coordinates
%                     (representing probability) of the PDF of the observed
%                     nearest neighbor distances between each subject event
%                     and their nearest landmark event. This must have as
%                     many elements as the input OBSERVED_PDF_X. Each
%                     element should have a value between 0 and 1, and the
%                     sum of all elements should equal 1.
%
% 3. OBSERVED_NEVENTS - Numeric scalar representing the sample size of
%                       subject events. This must be a positive integer.
%
% 4. RANDOM_PDF_X - Numeric column vector listing the x-coordinates
%                   (representing distance) of the PDF of the estimated
%                   random nearest neighbor distances between randomized
%                   subject events and their nearest landmark event. This
%                   must have as many elements as the intput RANDOM_PDF_Y
%                   and The elemenets must monotonically increase.
%
% 5. RANDOM_PDF_Y - Numeric column vector listing the x-coordinates
%                   (representing probability) of the PDF of the estimated
%                   random nearest neighbor distances between randomized
%                   subject events and their nearest landmark event. This
%                   must have as many elements as the input RANDOM_PDF_X.
%                   Each element should have a value between 0 and 1, and
%                   the sum of all elements should equal 1.
%
% 6. RANDOM_NEVENTS - Numeric scalar representing the sample size of
%                     randomized subject events. This must be a positive
%                     integer.
%
% 7. ALPHA - numeric scalar representing the significance level for
%            testing complete spatial randomness of the subject pattern
%            relative to the landmark pattern. This must be a  with a value
%            between 0 and 1, exlusive.
%
% -------------------------------------------------------------------------
%
% OUTPUTS: 
%
% 1. DELTA_CDF_X - Numeric column vector listing the x-coordinates
%                  (representing distance) of the delta function which
%                  compares the observed and random CDFs. The x-coordinates
%                  are calculated as the set of all unique x-coordinates
%                  from the observed and random PDFs. This will have as
%                  many elements as the output DELTA_CDF_Y and the
%                  elemenets will monotonically increase.
%
% 2. DELTA_CDF_Y - Numeric column vector listing the x-coordinates
%                  (representing probability) of the delta function which
%                  compares the observed and random CDFs. The y-coordinates
%                  are calculated as the random PDF y-coordinates
%                  subtracted from the observed PDF y-coordinates. This
%                  will have as many elements as the output DELTA_CDF_Y.
%                  Each element will have a value between -1 and 1.
%
% 3. TESTSTATISTIC - Numeric scalar that summarizes the magnitude of
%                    deviation from complete spatial randomness. Calculated
%                    as the value from DELTA_CDF_Y with the largest
%                    magnitude.
% 
% 4. PVALUE - Numeric scalar representing that asymptotic p-value
%             approximation from a 2-sided Kolmogorov-Smirnov test
%             comparing the observed and random CDFs. This wil have a value
%             between 0 and 1. If this has a value less than ALPHA, then
%             there is sufficient evidence to reject the nullhypothesis
%             which states that the observed and random distributions are
%             them same.
%
% 5. VERDICT - Logical scalar indicating whether the null hypothesis from
%              the Kolmogorov-Smirnov test was rejected. A value of 'true'
%              indicates there is sufficient evidence to reject the null
%              hypothesis and a value of 'false' indicates there is not
%              sufficient evidence to reject the null hypothesis.
%
% -------------------------------------------------------------------------

    if nargin < 7
        alpha = 0.05; % default value
    end
    
    % Convert to column vector to ensure consistent output format
    observed_PDF_x = observed_PDF_x(:);
    observed_PDF_y = observed_PDF_y(:);
    random_PDF_x = random_PDF_x(:);
    random_PDF_y = random_PDF_y(:);
    
    % Regenerate input distributions over a global x-coordinate scheme
    delta_CDF_x = unique([observed_PDF_x; random_PDF_x]); % find common bins
    [~, observedIndices] = ismember(observed_PDF_x, delta_CDF_x); % find where obs_pdf_y fits into global x-coordinates
    [~, randomIndices] = ismember(random_PDF_x, delta_CDF_x); % find where ran_pdf_y fits into global x-coordinates
    observed_PDF_y_global = zeros(size(delta_CDF_x)); % initialize new y-values as zeros
    random_PDF_y_global = observed_PDF_y_global; % initialize new y-values as zeros
    observed_PDF_y_global(observedIndices) = observed_PDF_y; % insert original y-values onto global x-coordinate scheme
    random_PDF_y_global(randomIndices) = random_PDF_y; % insert original y-values onto global x-coordinate scheme
    
    % Compute eCDFs from ePDFs
    observed_CDF_y = PDF2CDF(observed_PDF_y_global);
    random_CDF_y = PDF2CDF(random_PDF_y_global);
    
    % Compute the test statistic
    delta_CDF_y = observed_CDF_y - random_CDF_y;
    [ks_p, ks_idx] = max(abs(delta_CDF_y));
    testStatistic = delta_CDF_y(ks_idx);
    
    % Compute the asymptotic pValue approximation and accept or reject the
    % null hypothesis on the basis of the pValue (based on code found in
    % MATLAB's built-in kstest2 function)
    n = observed_nEvents * random_nEvents / (observed_nEvents + random_nEvents);
    lambda = max((sqrt(n) + 0.12 + 0.11 / sqrt(n)) * ks_p, 0);
    pValue = exp(-2 * lambda * lambda);
    verdict = (alpha >= pValue);

end

function [PDF_x, PDF_y] = ePDF(distances)
% Calculate the empirical probability density function (ePDF) of a set of
% distance measurements. Optionally, events can be weighted, if for example
% each event occupies a different volume of space.
%
% 
% -------------------------------------------------------------------------
%
% INPUTS:
%
% 1. DISTANCES - Numeric column vector containing samples of distance
%                measurements which will be used to generate a PDF.
%
% -------------------------------------------------------------------------
%
% OUTPUTS: 
%
% 1. PDF_X - Numeric column vector of monotonically increasing values
%            representing the x-coordinates (bin values) of the PDF which
%            describes the distribution of the input sample distance
%            measurements. This is equavalent to all of the unique values
%            from DISTANCES sorted in ascending order. This will have the
%            same number of elements as the output PDF_Y.
%
% 2. PDF_Y - Numeric column vector of y-coordinates (bin counts) of the PDF
%            which describes the distribution of the input sample distance
%            measurements. This is equavalent to the number of distance
%            samples which have the corresponding value in PDF_X. This will
%            have the same number of elements as the output PDF_X.
%
% -------------------------------------------------------------------------

    [PDF_y, PDF_x] = groupcounts(distances);
    PDF_x = double(PDF_x);

    % make sure PDF is defined for full distance range [0, maxDistance]
    if PDF_x(1) ~= 0
        PDF_x = [0; PDF_x];
        PDF_y = [0; PDF_y];
    end

end

function CDF_y = PDF2CDF(PDF_y)
% Compute an empirical cumulative distribution function (eCDF) based on
% an input empirical probability density function (ePDF).
%
% 
% -------------------------------------------------------------------------
%
% INPUTS:
%
% 1. PDF_Y - Numeric column vector of y-coordinates (bin counts) of a PDF
%            which describes the distribution of some sampling of
%            measurements. This will have the same number of elements as
%            the output CDF_Y.
%
%
% -------------------------------------------------------------------------
%
% OUTPUTS: 
%
% 1. CDF_Y - Numeric column vector of monotonically increasing
%            y-coordinates (cumulative probabilties) of a CDF which
%            describes the distribution of some sampling of measurements.
%            This will have the same number of elements as the input PDF_Y
%            and the values will be between 0 and 1 and the last element is
%            garaunteed to be 1.
%
%
% -------------------------------------------------------------------------

    CDF_y = cumsum(PDF_y);
    CDF_y = CDF_y ./ CDF_y(end);

end
