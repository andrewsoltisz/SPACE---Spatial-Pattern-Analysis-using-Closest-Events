# SPACE---Spatial-Pattern-Analysis-using-Closest-Events
Statistically characterize the spatial relationship between 2 image components using point pattern analysis.  

## Authorship:  
**Title**: SPACE - Spatial Pattern Analysis using Closest Events  
**Citation**: Andrew M Soltisz, Peter F Craigmile, Rengasayee Veeraraghavan. Spatial Pattern Analysis using Closest Events (SPACE)—A Nearest Neighbor Point Pattern Analysis Framework for Assessing Spatial Relationships from Digital Images. Microscopy and Microanalysis, 2024, ozae022. https://doi.org/10.1093/mam/ozae022.  
**Author**: Andrew M. Soltisz  
**Contact**: andysoltisz@gmail.com  
**LinkedIn**: https://www.linkedin.com/in/andrew-soltisz/  
**Acknowledgements**: Rengasayee Veeraraghavan, Seth H. Weinberg, Peter F. Craigmile  
**Theory**: SPACE leverages bivariate nearest neighor distribution functions (G- and F-functions) from the point pattern analysis paradigm. The underlying theory is described in https://doi.org/10.1002/9780470057339.van007.  

Copyright (C) 2025, Andrew Michael Soltisz.  
All rights reserved.   
The source code within this directory is licensed under the BSD-3-Clause License found in the LICENSE.txt file.

 
## Description:  
Spatial Pattern Analysis using Closest Events (SPACE) is an image analysis framework for characterizing the spatial relationship between 2 image components using nearest neighbor-based spatial statistics. Given 2 image components X and Y, identified by the true pixels in their respective image masks X_mask and Y_mask, the spatial relationship of X relative to Y (X-->Y) and Y relative to X (Y-->X) are characterized by comparing the distribution of inter-pixel nearest neighbor distances (“observed” G-function) to the distribution produced by the distances between every pixel in the image to the nearest component pixel (“random” F-function). This relationship is summarized using a single output parameter, the spatial association index, computed as the maximum absolute difference between these functions, where negative values indicate a dispersed spatial relationship, positive values indicate aggregation, and values at-or-near 0 indicate that the image components are completely spatially random relative to one another. The function output will be a table whose fields specify the coordinates of the aforementioned empirical distribution functions and resulting spatial association index. Batch analysis of a set of SPACE results can be performed to summarize the global statistical behaviour across multiple image pairs. The output will be a table whose fields specify the coordinates of the median-versions of the aforementioned empirical distribution functions and resulting global spatial association index. Additionally, you can specify the analysis be performed on only a subregion of the image by including a third input indicating the region of interest as a logical mask. Example scripts and images are provided to demonstrate the use of this framework.


## Included Files:  
SPACE.m - the main function for performing SPACE analysis.  
SPACE_batch.m - function for performing batch analysis of SPACE results.  
ttest2w.m - accessory function to compute the weighted version of a 2-sided Student's T-test. Recommended test for comparing spatial association indices across groups of images to account for differences in signal volume.  
plot_error.m - accessory function to plot lines with shaded error region.  
gen_overlay.m - function to combined masks into a single RGB image.  
example_single.m - script with an example use of SPACE to analyze a single image.  
example_batch.m - script with an example use of SPACE to analyze a batch of images.  
example_compare_groups.m - script with an example use of SPACE to analyze and compare 3 different groups of images.  
generate_example_data.m - script called by example scripts to generate data for demonstration of SPACE.  

## Abbreviations:  
CDF – empirical cumulative distribution function  
PDF – empirical probability density function  
NN – nearest neighbor  
SPACE – Spatial Pattern Analysis using Closest Events  
KS – nonparametric 2-sided Kolmogorov-Smirnov test of equality between 2 empirical distributions  
CSR – complete spatial randomness  
ROI – region of interest  
 
## Table field definitions for SPACE single-image analysis function output. Rows correspond to individual images, columns correspond to SPACE results for each image:  
1. **Observed_Event_Count** – number of events (white pixels) used to generate the observed CDF.  
2. **Random_Event_Count** – number of points (pixels) used to generate the random CDF.  
3. **Observed_x** – x-coordinates (distance) for all observed functions.  
4. **Observed_PDF_y** – y-coordinates (event count) for PDF of observed event NN distances. Corresponding x-coordinates are stored in Observed_x.  
5. **Random_x** – x-coordinates (distance) for all random functions.  
6. **Random_PDF_y** – y-coordinates (event count) for PDF of random event distances. Corresponding x-coordinates are stored in Random_x.  
7. **Observed_CDF_y** – y-coordinates (probability) for CDF of observed event NN distances. Corresponding x-coordinates are stored in Observed_x.  
8. **Random_CDF_y** – y-coordinates (probability) for CDF of random event NN distances. Corresponding x-coordinates are stored in Random_x.  
9. **Delta_CDF_x** – x-coordinates (distance) of delta function. Delta functions represent subtraction of the random CDF from the observed CDF.  
10. **Delta_CDF_y** – y-coordinates (delta probability) of delta function. Corresponding x-coordinates are stored in Delta_CDF_x.  
11. **Spatial_Association_Index** – spatial association index for describing that type and magnitude of spatial association exhibited between both signals. Computed as the absolute maximum of XY_Delta_CDF_y. This can take on a value between -1 and 1, where negative values indicate dispersion of signals between both masks, positive values indicate aggregation, and values at or near zero indicate random (Poisson) association.  
12. **Spatial_Association_pValue** – p-value from a 2-sided KS test comparing the observed and random distributions. The null hypothesis states that these distributions are identical and therefor signals from the first mask are CSR relative to signals from the second mask.
13. **Spatial_Association_Verdict** - logical scalar indicating whether the null hypothesis from the KS test was rejected. A value of 'true' indicates there is sufficient evidence to reject the null hypothesis and a value of 'false' indicates there is not sufficient evidence to reject the null hypothesis.
 
## Table field definitions for SPACE batch analysis function output. Rows correspond to individual sets of SPACE analyses, columns correspond to diferent results from the batch analysis:  
1. **Sample_Size** – number of images used to generate the SPACE results.  
2. **Global_x** – x-coordinates shared by all functions and used to create their global y-coordinates.  
3. **Observed_CDF_x** – subset of global x-coordinates (Global_x) corresponding to observed CDF y-coordinates (Observed_CDF_y_Matrix). This is a single column vector used for all observed CDFs.  
4. **Observed_CDF_y_Matrix** – y-coordinates for all individual observed CDFs that are all defined over the same global x-coordinate scheme (Observed_CDF_x). Functions for individual images are stored in every column and y-coordinates are stored along rows. Thus, the number of columns corresponds to the number of individual analyses and the number of rows corresponds to the number of y-coordinates for each function.  
5. **Random_CDF_x** – subset of global x-coordinates (Global_x) corresponding to random CDF y-coordinates (Random_CDF_y_Matrix). This is a single column vector used for all random CDFs.  
6. **Random_CDF_y_Matrix** – y-coordinates for all individual random CDFs that are all defined over the same global x-coordinate scheme (Random_CDF_x). Functions for individual analyses are stored in every column and y-coordinates are stored along rows. Thus, the number of columns corresponds to the number of individual analyses and the number of rows corresponds to the number of y-coordinates for each function.  
7. **Delta_CDF_x** – subset of global x-coordinates (Global_x) corresponding to delta function y-coordinates (Delta_CDF_y_Matrix). This is a single column vector used for all delta functions. Delta functions represent subtraction of the random CDF from the observed CDF.  
8. **Delta_CDF_y_Matrix** – y-coordinates for all individual delta functions that are all defined over the same global x-coordinate scheme (Delta_CDF_x). Functions for individual analyses are stored in every column and y-coordinates are stored along rows. Thus, the number of columns corresponds to the number of individual analyses and the number of rows corresponds to the number of y-coordinates for each function.  
9. **Observed_CDF_y_Median** – y-coordinates for the median observed CDF. Computed as the weighted median for every row of Observed_CDF_y_Matrix. Corresponding x-coordinates are stored in Observed_CDF_x.   
10. **Observed_CDF_y_Lower** – y-coordinates for the lower quantile of the median observed CDF. Computed as the weighted lower quantile for every row of Observed_CDF_y_Matrix. Corresponding x-coordinates are stored in Observed_CDF_x.  
11. **Observed_CDF_y_Upper** – y-coordinates for the upper quantile of the median observed CDF. Computed as the weighted upper quantile for every row of Observed_CDF_y_Matrix. Corresponding x-coordinates are stored in Observed_CDF_x.  
12. **Random_CDF_y_Median** – y-coordinates for the median random CDF. Computed as the weighted median for every row of Random_CDF_y_Matrix. Corresponding x-coordinates are stored in Random_CDF_x.  
13. **Random_CDF_y_Lower** – y-coordinates for the lower quantile of the median random CDF. Computed as the weighted lower quantile for every row of Random_CDF_y_Matrix. Corresponding x-coordinates are stored in Random_CDF_x.  
14. **Random_CDF_y_Upper** – y-coordinates for the upper quantile of the median random CDF. Computed as the weighted upper quantile for every row of Random_CDF_y_Matrix. Corresponding x-coordinates are stored in Random_CDF_x.  
15. **Delta_CDF_y_Median** – y-coordinates for the median delta function. Computed as the weighted median for every row of Delta_CDF_y_Matrix. Corresponding x-coordinates are stored in Delta_CDF_x.  
16. **Delta_CDF_y_Lower** – y-coordinates for the lower quantile of the median delta function. Computed as the weighted lower quantile for every row of Delta_CDF_y_Matrix. Corresponding x-coordinates are stored in Delta_CDF_x.  
17. **Delta_CDF_y_Upper** – y-coordinates for the upper quantile of the median delta function. Computed as the weighted upper quantile for every row of Delta_CDF_y_Matrix. Corresponding x-coordinates are stored in Delta_CDF_x.  
18. **Sample_Event_Count** – List of event counts (white pixels) for each individual image, used to generate all functions. This is a column vector and each element corresponds to the columns in all function matrices.   
19. **Sample_Spatial_Association_Index** – List of spatial association indices from individual analyses. Computed as the absolute maximum of each column in Delta_CDF_y_Matrix.   
20. **Global_Spatial_Association_Index** – Global spatial association index given all analyses in the batch. Computed as the absolute maximum of Delta_CDF_y_Median.  
21. **Global_CSR_Verdict** – verdict of the test for CSR between the first and second mask. The null hypothesis states that the white pixels in one mask are CSR relative to the white pixels in the other mask. This variable returns true if the null hypothesis is rejected and false if it is not rejected.
