# SPACE---Spatial-Pattern-Analysis-using-Closest-Events
Statistically characterize the spatial relationship between 2 image components using point pattern analysis.

Title: SPACE - Spatial Pattern Analysis using Closest Events
Author: Andrew M. Soltisz
Contact: andysoltisz@gmail.com
LinkedIn: https://www.linkedin.com/in/andrew-soltisz/
Acknowledgements: Rengasayee Veeraraghavan, Seth H. Weinberg, Peter F. Craigmile
 
Description: 
Spatial Pattern Analysis using Closest Events (SPACE) is an image analysis framework for characterizing the spatial relationship between 2 image components using nearest neighbor-based spatial statistics. Given 2 image components X and Y, identified by the true pixels in their respective image masks X_mask and Y_mask, the spatial relationship of X relative to Y (X-->Y) and Y relative to X (Y-->X) are characterized by comparing the distribution of inter-pixel nearest neighbor distances (“observed” G-function) to the distribution produced by the distances between every pixel in the image to the nearest component pixel (“random” F-function). This relationship is summarized using a single output parameter, the spatial association index, computed as the maximum absolute difference between these functions, where negative values indicate a dispersed spatial relationship, positive values indicate aggregation, and values at-or-near 0 indicate that the image components are completely spatially random relative to one another. A single image can be analyzed by directly providing both logical masks as separate inputs to the SPACE function. The output will be a table whose fields specify the coordinates of the aforementioned empirical distribution functions and resulting spatial association index.  Batch analysis can be performed on a group of images by inputting cell arrays of logical masks. The output will be a table whose fields specify the coordinates of the median-versions of the aforementioned empirical distribution functions and resulting global spatial association index.  Additionally, you can specify the analysis be performed on only a subregion of the image by including a third input indicating the region of interest using as a logical mask. Example scripts and images are provided to demonstrate the use of this framework.

Included Files:
SPACE.m - the main function for performing SPACE analysis. 
ttest2w.m - accessory function to compute the weighted version of a 2-sided Student's T-test. Recommended test for comparing spatial association indices across groups of images to account for differences in signal volume.
plot_error.m - accessory function to plot lines with shaded error region.
example_single.m - script with an example use of SPACE to analyze a single image.
example_batch.m - script with an example use of SPACE to analyze a batch of images.
example_compare_groups.m - script with an example use of SPACE to analyze and compare 3 different groups of images.
generate_example_data.m - script called by example scripts to generate data for demonstration of SPACE.
Abbreviation Definitions:
CDF – empirical cumulative distribution function
PDF – empirical probability density function
NN – nearest neighbor
SPACE – Spatial Pattern Analysis using Closest Events
KS – nonparametric 2-sided Kolmogorov-Smirnov test of equality between 2 empirical distributions
CSR – complete spatial randomness
ROI – region of interest
 
Table field definitions for SPACE single-image analysis function output. Rows correspond to individual images, columns correspond to SPACE results for each image:
XY_Observed_Event_Count – number of X-events used to generate the X-->Y observed CDF
XY_Random_Event_Count – number of points used to generate the X-->Y random CDF
XY_Observed_x – x-coordinates (distance) for all observed X-->Y functions
XY_Observed_PDF_y – y-coordinates (event count) for PDF of observed X-->Y event NN distances
XY_Random_x – x-coordinates (distance) for all random X-->Y functions
XY_Random_PDF_y – y-coordinates (event count) for PDF of random X-->Y event distances
XY_Observed_CDF_y – y-coordinates (probability) for CDF of observed X-->Y event NN distances
XY_Random_CDF_y – y-coordinates (probability) for CDF of random X-->Y event NN distances
XY_Delta_CDF_x – x-coordinates (distance) of X-->Y delta function
XY_Delta_CDF_y – y-coordinates (delta probability) of X-->Y delta function
XY_Spatial_Association_Index – spatial association index for X-->Y. Computed as the absolute maximum of XY_Delta_CDF_y.
XY_Spatial_Assocation_pValue – p-value from a 2-sided KS test comparing the X-->Y observed and random distributions. The null hypothesis states that these distributions are identical and therefor X-->Y is CSR
YX_Observed_Event_Count – number of Y-events used to generate the Y-->X observed CDF
YX_Random_Event_Count – number of points used to generate the Y-->X random CDF
YX_Observed_x – x-coordinates (distance) for all observed Y-->X functions
YX_Observed_PDF_y – y-coordinates (event count) for PDF of observed Y-->X event NN distances
YX_Random_x – x-coordinates (distance) for all random Y-->X functions
YX_Random_PDF_y – y-coordinates (event count) for PDF of random Y-->X event distances
YX_Observed_CDF_y – y-coordinates (probability) for CDF of observed Y-->X event NN distances
YX_Random_CDF_y – y-coordinates (probability) for CDF of random Y-->X event NN distances
YX_Delta_CDF_x – x-coordinates (distance) of Y-->X delta function
YX_Delta_CDF_y – y-coordinates (delta probability) of Y-->X delta function
YX_Spatial_Association_Index – spatial association index for Y-->X. Computed as the absolute maximum of YX_Delta_CDF_y. 
YX_Spatial_Assocation_pValue – p-value from a 2-sided KS test comparing the Y-->X observed and random distributions. The null hypothesis states that these distributions are identical and therefor X-->Y is CSR
 
Table field definitions for SPACE batch-image analysis function output. Rows correspond to individual batches of images, columns correspond to batch SPACE results for each group of images:
Sample_Size – number of images in the batch
XY_Global_x – x-coordinates shared by all X-->Y functions and used to create their global y-coordinates
XY_Observed_CDF_x – subset of global x-coordinates (XY_Global_x) corresponding to observed X-->Y CDF y-coordinates (XY_Observed_CDF_y_Matrix). This is a single column vector used for all observed X-->Y CDFs.
XY_Observed_CDF_y_Matrix – y-coordinates for all individual observed X-->Y CDFs that are all defined over the same global x-coordinate scheme (XY_Observed_CDF_x). Functions for individual images are stored in every column and y-coordinates are stored along rows. Thus, the number of columns corresponds to the number of individual images and the number of rows corresponds to the number of y-coordinates for each function.
XY_Random_CDF_x – subset of global x-coordinates (XY_Global_x) corresponding to random X-->Y CDF y-coordinates (XY_Random_CDF_y_Matrix). This is a single column vector used for all random X-->Y CDFs.
XY_Random_CDF_y_Matrix – y-coordinates for all individual random X-->Y CDFs that are all defined over the same global x-coordinate scheme (XY_Random_CDF_x). Functions for individual images are stored in every column and y-coordinates are stored along rows. Thus, the number of columns corresponds to the number of individual images and the number of rows corresponds to the number of y-coordinates for each function.
XY_Delta_CDF_x – subset of global x-coordinates (XY_Global_x) corresponding to X-->Y delta function y-coordinates (XY_Delta_CDF_y_Matrix). This is a single column vector used for all X-->Y delta functions.
XY_Delta_CDF_y_Matrix – y-coordinates for all individual X-->Y delta functions that are all defined over the same global x-coordinate scheme (XY_Delta_CDF_x). Functions for individual images are stored in every column and y-coordinates are stored along rows. Thus, the number of columns corresponds to the number of individual images and the number of rows corresponds to the number of y-coordinates for each function.
XY_Observed_CDF_y_Median – y-coordinates for the median X-->Y observed CDF. Computed as the weighted median for every row of XY_Observed_CDF_y_Matrix. Corresponding x-coordinates are stored in XY_Observed_CDF_x. 
XY_Observed_CDF_y_Lower – y-coordinates for the lower quantile of the median X-->Y observed CDF. Computed as the weighted lower quantile for every row of XY_Observed_CDF_y_Matrix. Corresponding x-coordinates are stored in XY_Observed_CDF_x.
XY_Observed_CDF_y_Upper – y-coordinates for the upper quantile of the median X-->Y observed CDF. Computed as the weighted upper quantile for every row of XY_Observed_CDF_y_Matrix. Corresponding x-coordinates are stored in XY_Observed_CDF_x.
XY_Random_CDF_y_Median – y-coordinates for the median X-->Y random CDF. Computed as the weighted median for every row of XY_Random_CDF_y_Matrix. Corresponding x-coordinates are stored in XY_Random_CDF_x.
XY_Random_CDF_y_Lower – y-coordinates for the lower quantile of the median X-->Y random CDF. Computed as the weighted lower quantile for every row of XY_Random_CDF_y_Matrix. Corresponding x-coordinates are stored in XY_Random_CDF_x.
XY_Random_CDF_y_Upper – y-coordinates for the upper quantile of the median X-->Y random CDF. Computed as the weighted upper quantile for every row of XY_Random_CDF_y_Matrix. Corresponding x-coordinates are stored in XY_Random_CDF_x.
XY_Delta_CDF_y_Median – y-coordinates for the median X-->Y delta function. Computed as the weighted median for every row of XY_Delta_CDF_y_Matrix. Corresponding x-coordinates are stored in XY_Delta_CDF_x.
XY_Delta_CDF_y_Lower – y-coordinates for the lower quantile of the median X-->Y delta function. Computed as the weighted lower quantile for every row of XY_Delta_CDF_y_Matrix. Corresponding x-coordinates are stored in XY_Delta_CDF_x.
XY_Delta_CDF_y_Upper – y-coordinates for the upper quantile of the median X-->Y delta function. Computed as the weighted upper quantile for every row of XY_Delta_CDF_y_Matrix. Corresponding x-coordinates are stored in XY_Delta_CDF_x.
XY_Sample_Event_Count – List of X-event counts for each individual image, used to generate all X-->Y functions. This is a column vector and each element corresponds to the columns in all X-->Y function matrices. 
XY_Sample_Spatial_Association_Index – List of X-->Y spatial association indices for individual images. Computed as the absolute maximum of each column in XY_Delta_CDF_y_Matrix. 
XY_Global_Spatial_Association_Index – Global spatial association index of X-->Y given all images in the batch. Computed as the absolute maximum of XY_Delta_CDF_y_Median.
XY_Global_CSR_Verdict – verdict of the test for CSR in X-->Y. Equal to ‘true’ if the X-->Y median delta function quantile envelope (XY_Delta_CDF_y_Lower and XY_Delta_CDF_y_Upper) ever deviates from 0, and equal to ‘false’ if the envelope overlaps with 0 for the full distance range.
YX_Global_x – x-coordinates shared by all Y-->X functions and used to create their global y-coordinates
YX_Observed_CDF_x – subset of global x-coordinates (YX_Global_x) corresponding to observed Y-->X CDF y-coordinates (YX_Observed_CDF_y_Matrix). This is a single column vector used for all observed Y-->X CDFs. 
YX_Observed_CDF_y_Matrix – y-coordinates for all individual observed Y-->X CDFs that are all defined over the same global x-coordinate scheme (YX_Observed_CDF_x). Functions for individual images are stored in every column and y-coordinates are stored along rows. Thus, the number of columns corresponds to the number of individual images and the number of rows corresponds to the number of y-coordinates for each function.
YX_Random_CDF_x – subset of global x-coordinates (YX_Global_x) corresponding to random Y-->X CDF y-coordinates (YX_Random_CDF_y_Matrix). This is a single column vector used for all random YàX CDFs. 
YX_Random_CDF_y_Matrix – y-coordinates for all individual random Y-->X CDFs that are all defined over the same global x-coordinate scheme (YX_Random_CDF_x). Functions for individual images are stored in every column and y-coordinates are stored along rows. Thus, the number of columns corresponds to the number of individual images and the number of rows corresponds to the number of y-coordinates for each function.
YX_Delta_CDF_x – subset of global x-coordinates (YX_Global_x) corresponding to Y-->X delta function y-coordinates (YX_Delta_CDF_y_Matrix). This is a single column vector used for all Y-->X delta functions. 
YX_Delta_CDF_y_Matrix – y-coordinates for all individual Y-->X delta functions that are all defined over the same global x-coordinate scheme (YX_Delta_CDF_x). Functions for individual images are stored in every column and y-coordinates are stored along rows. Thus, the number of columns corresponds to the number of individual images and the number of rows corresponds to the number of y-coordinates for each function.
YX_Observed_CDF_y_Median – y-coordinates for the median Y-->X observed CDF. Computed as the weighted median for every row of YX_Observed_CDF_y_Matrix. Corresponding x-coordinates are stored in YX_Observed_CDF_x. 
YX_Observed_CDF_y_Lower – y-coordinates for the lower quantile of the median Y-->X observed CDF. Computed as the weighted lower quantile for every row of YX_Observed_CDF_y_Matrix. Corresponding x-coordinates are stored in YX_Observed_CDF_x.
YX_Observed_CDF_y_Upper – y-coordinates for the upper quantile of the median Y-->X observed CDF. Computed as the weighted upper quantile for every row of YX_Observed_CDF_y_Matrix. Corresponding x-coordinates are stored in YX_Observed_CDF_x.
YX_Random_CDF_y_Median – y-coordinates for the median Y-->X random CDF. Computed as the weighted median for every row of YX_Random_CDF_y_Matrix. Corresponding x-coordinates are stored in YX_Random_CDF_x.
YX_Random_CDF_y_Lower – y-coordinates for the lower quantile of the median Y-->X random CDF. Computed as the weighted lower quantile for every row of YX_Random_CDF_y_Matrix. Corresponding x-coordinates are stored in YX_Random_CDF_x.
YX_Random_CDF_y_Upper – y-coordinates for the upper quantile of the median Y-->X random CDF. Computed as the weighted upper quantile for every row of YX_Random_CDF_y_Matrix. Corresponding x-coordinates are stored in YX_Random_CDF_x.
YX_Delta_CDF_y_Median – y-coordinates for the median Y-->X delta function. Computed as the weighted median for every row of YX_Delta_CDF_y_Matrix. Corresponding x-coordinates are stored in YX_Delta_CDF_x.
YX_Delta_CDF_y_Lower – y-coordinates for the lower quantile of the median Y-->X delta function. Computed as the weighted lower quantile for every row of YX_Delta_CDF_y_Matrix. Corresponding x-coordinates are stored in YX_Delta_CDF_x.
YX_Delta_CDF_y_Upper – y-coordinates for the upper quantile of the median Y-->X delta function. Computed as the weighted upper quantile for every row of YX_Delta_CDF_y_Matrix. Corresponding x-coordinates are stored in YX_Delta_CDF_x.
YX_Sample_Event_Count – List of X-event counts for each individual image, used to generate all Y-->X functions. This is a column vector and each element corresponds to the columns in all Y-->X function matrices. 
YX_Sample_Spatial_Association_Index – List of Y-->X spatial association indices for individual images. Computed as the absolute maximum of each column in YX_Delta_CDF_y_Matrix. 
YX_Global_Spatial_Association_Index – Global spatial association index of Y-->X given all images in the batch. Computed as the absolute maximum of YX_Delta_CDF_y_Median.
YX_Global_CSR_Verdict – verdict of the test for CSR in Y-->X. Equal to ‘true’ if the Y-->X median delta function quantile envelope (YX_Delta_CDF_y_Lower and YX_Delta_CDF_y_Upper) ever deviates from 0, and equal to ‘false’ if the envelope overlaps with 0 for the full distance range.
