# SPACE---Spatial-Pattern-Analysis-using-Closest-Events
Statistically characterize the spatial relationship between 2 image components using point pattern analysis.  

## Authorship:  
**Title**: SPACE - Spatial Pattern Analysis using Closest Events  
**Author**: Andrew M. Soltisz  
**Contact**: andysoltisz@gmail.com  
**LinkedIn**: https://www.linkedin.com/in/andrew-soltisz/  
**Acknowledgements**: Rengasayee Veeraraghavan, Seth H. Weinberg, Peter F. Craigmile  
**Theory**: SPACE leverages bivariate nearest neighor distribution functions (G- and F-functions) from the point pattern analysis paradigm. The underlying theory is described in https://doi.org/10.1002/9780470057339.van007.  

 
## Description:  
Spatial Pattern Analysis using Closest Events (SPACE) is an image analysis framework for characterizing the spatial relationship between 2 image components using nearest neighbor-based spatial statistics. Given 2 image components X and Y, identified by the true pixels in their respective image masks X_mask and Y_mask, the spatial relationship of X relative to Y (X-->Y) and Y relative to X (Y-->X) are characterized by comparing the distribution of inter-pixel nearest neighbor distances (“observed” G-function) to the distribution produced by the distances between every pixel in the image to the nearest component pixel (“random” F-function). This relationship is summarized using a single output parameter, the spatial association index, computed as the maximum absolute difference between these functions, where negative values indicate a dispersed spatial relationship, positive values indicate aggregation, and values at-or-near 0 indicate that the image components are completely spatially random relative to one another. A single image can be analyzed by directly providing both logical masks as separate inputs to the SPACE function. The output will be a table whose fields specify the coordinates of the aforementioned empirical distribution functions and resulting spatial association index.  Batch analysis can be performed on a group of images by inputting cell arrays of logical masks. The output will be a table whose fields specify the coordinates of the median-versions of the aforementioned empirical distribution functions and resulting global spatial association index.  Additionally, you can specify the analysis be performed on only a subregion of the image by including a third input indicating the region of interest using as a logical mask. Example scripts and images are provided to demonstrate the use of this framework.


## Included Files:  
SPACE.m - the main function for performing SPACE analysis.  
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
1. **XY_Observed_Event_Count** – number of X-events used to generate the X-->Y observed CDF  
2. **XY_Random_Event_Count** – number of points used to generate the X-->Y random CDF  
3. **XY_Observed_x** – x-coordinates (distance) for all observed X-->Y functions  
4. **XY_Observed_PDF_y** – y-coordinates (event count) for PDF of observed X-->Y event NN distances  
5. **XY_Random_x** – x-coordinates (distance) for all random X-->Y functions  
6. **XY_Random_PDF_y** – y-coordinates (event count) for PDF of random X-->Y event distances  
7. **XY_Observed_CDF_y** – y-coordinates (probability) for CDF of observed X-->Y event NN distances  
8. **XY_Random_CDF_y** – y-coordinates (probability) for CDF of random X-->Y event NN distances  
9. **XY_Delta_CDF_x** – x-coordinates (distance) of X-->Y delta function  
10. **XY_Delta_CDF_y** – y-coordinates (delta probability) of X-->Y delta function  
11. **XY_Spatial_Association_Index** – spatial association index for X-->Y. Computed as the absolute maximum of XY_Delta_CDF_y.  
12. **XY_Spatial_Assocation_pValue** – p-value from a 2-sided KS test comparing the X-->Y observed and random distributions. The null hypothesis states that these distributions are identical and therefor X-->Y is CSR  
13. **YX_Observed_Event_Count** – number of Y-events used to generate the Y-->X observed CDF  
14. **YX_Random_Event_Count** – number of points used to generate the Y-->X random CDF  
15. **YX_Observed_x** – x-coordinates (distance) for all observed Y-->X functions  
16. **YX_Observed_PDF_y** – y-coordinates (event count) for PDF of observed Y-->X event NN distances  
17. **YX_Random_x** – x-coordinates (distance) for all random Y-->X functions  
18. **YX_Random_PDF_y** – y-coordinates (event count) for PDF of random Y-->X event distances  
19. **YX_Observed_CDF_y** – y-coordinates (probability) for CDF of observed Y-->X event NN distances  
20. **YX_Random_CDF_y** – y-coordinates (probability) for CDF of random Y-->X event NN distances  
21. **YX_Delta_CDF_x** – x-coordinates (distance) of Y-->X delta function  
22. **YX_Delta_CDF_y** – y-coordinates (delta probability) of Y-->X delta function  
23. **YX_Spatial_Association_Index** – spatial association index for Y-->X. Computed as the absolute maximum of YX_Delta_CDF_y.   
24. **YX_Spatial_Assocation_pValue** – p-value from a 2-sided KS test comparing the Y-->X observed and random distributions. The null hypothesis states that these distributions are identical and therefor X-->Y is CSR  
 
## Table field definitions for SPACE batch-image analysis function output. Rows correspond to individual batches of images, columns correspond to batch SPACE results for each group of images:  
1. **Sample_Size** – number of images in the batch  
2. **XY_Global_x** – x-coordinates shared by all X-->Y functions and used to create their global y-coordinates  
3. **XY_Observed_CDF_x** – subset of global x-coordinates (XY_Global_x) corresponding to observed X-->Y CDF y-coordinates (XY_Observed_CDF_y_Matrix). This is a single column vector used for all observed X-->Y CDFs.  
4. **XY_Observed_CDF_y_Matrix** – y-coordinates for all individual observed X-->Y CDFs that are all defined over the same global x-coordinate scheme (XY_Observed_CDF_x). Functions for individual images are stored in every column and y-coordinates are stored along rows. Thus, the number of columns corresponds to the number of individual images and the number of rows corresponds to the number of y-coordinates for each function.  
5. **XY_Random_CDF_x** – subset of global x-coordinates (XY_Global_x) corresponding to random X-->Y CDF y-coordinates (XY_Random_CDF_y_Matrix). This is a single column vector used for all random X-->Y CDFs.  
6. **XY_Random_CDF_y_Matrix** – y-coordinates for all individual random X-->Y CDFs that are all defined over the same global x-coordinate scheme (XY_Random_CDF_x). Functions for individual images are stored in every column and y-coordinates are stored along rows. Thus, the number of columns corresponds to the number of individual images and the number of rows corresponds to the number of y-coordinates for each function.  
7. **XY_Delta_CDF_x** – subset of global x-coordinates (XY_Global_x) corresponding to X-->Y delta function y-coordinates (XY_Delta_CDF_y_Matrix). This is a single column vector used for all X-->Y delta functions.  
8. **XY_Delta_CDF_y_Matrix** – y-coordinates for all individual X-->Y delta functions that are all defined over the same global x-coordinate scheme (XY_Delta_CDF_x). Functions for individual images are stored in every column and y-coordinates are stored along rows. Thus, the number of columns corresponds to the number of individual images and the number of rows corresponds to the number of y-coordinates for each function.  
9. **XY_Observed_CDF_y_Median** – y-coordinates for the median X-->Y observed CDF. Computed as the weighted median for every row of XY_Observed_CDF_y_Matrix. Corresponding x-coordinates are stored in XY_Observed_CDF_x.   
10. **XY_Observed_CDF_y_Lower** – y-coordinates for the lower quantile of the median X-->Y observed CDF. Computed as the weighted lower quantile for every row of XY_Observed_CDF_y_Matrix. Corresponding x-coordinates are stored in XY_Observed_CDF_x.  
11. **XY_Observed_CDF_y_Upper** – y-coordinates for the upper quantile of the median X-->Y observed CDF. Computed as the weighted upper quantile for every row of XY_Observed_CDF_y_Matrix. Corresponding x-coordinates are stored in XY_Observed_CDF_x.  
12. **XY_Random_CDF_y_Median** – y-coordinates for the median X-->Y random CDF. Computed as the weighted median for every row of XY_Random_CDF_y_Matrix. Corresponding x-coordinates are stored in XY_Random_CDF_x.  
13. **XY_Random_CDF_y_Lower** – y-coordinates for the lower quantile of the median X-->Y random CDF. Computed as the weighted lower quantile for every row of XY_Random_CDF_y_Matrix. Corresponding x-coordinates are stored in XY_Random_CDF_x.  
14. **XY_Random_CDF_y_Upper** – y-coordinates for the upper quantile of the median X-->Y random CDF. Computed as the weighted upper quantile for every row of XY_Random_CDF_y_Matrix. Corresponding x-coordinates are stored in XY_Random_CDF_x.  
15. **XY_Delta_CDF_y_Median** – y-coordinates for the median X-->Y delta function. Computed as the weighted median for every row of XY_Delta_CDF_y_Matrix. Corresponding x-coordinates are stored in XY_Delta_CDF_x.  
16. **XY_Delta_CDF_y_Lower** – y-coordinates for the lower quantile of the median X-->Y delta function. Computed as the weighted lower quantile for every row of XY_Delta_CDF_y_Matrix. Corresponding x-coordinates are stored in XY_Delta_CDF_x.  
17. **XY_Delta_CDF_y_Upper** – y-coordinates for the upper quantile of the median X-->Y delta function. Computed as the weighted upper quantile for every row of XY_Delta_CDF_y_Matrix. Corresponding x-coordinates are stored in XY_Delta_CDF_x.  
18. **XY_Sample_Event_Count** – List of X-event counts for each individual image, used to generate all X-->Y functions. This is a column vector and each element corresponds to the columns in all X-->Y function matrices.   
19. **XY_Sample_Spatial_Association_Index** – List of X-->Y spatial association indices for individual images. Computed as the absolute maximum of each column in XY_Delta_CDF_y_Matrix.   
20. **XY_Global_Spatial_Association_Index** – Global spatial association index of X-->Y given all images in the batch. Computed as the absolute maximum of XY_Delta_CDF_y_Median.  
21. **XY_Global_CSR_Verdict** – verdict of the test for CSR in X-->Y. Equal to ‘true’ if the X-->Y median delta function quantile envelope (XY_Delta_CDF_y_Lower and XY_Delta_CDF_y_Upper) ever deviates from 0, and equal to ‘false’ if the envelope overlaps with 0 for the full distance range.  
22. **YX_Global_x** – x-coordinates shared by all Y-->X functions and used to create their global y-coordinates  
23. **YX_Observed_CDF_x** – subset of global x-coordinates (YX_Global_x) corresponding to observed Y-->X CDF y-coordinates (YX_Observed_CDF_y_Matrix). This is a single column vector used for all observed Y-->X CDFs.   
24. **YX_Observed_CDF_y_Matrix** – y-coordinates for all individual observed Y-->X CDFs that are all defined over the same global x-coordinate scheme (YX_Observed_CDF_x). Functions for individual images are stored in every column and y-coordinates are stored along rows. Thus, the number of columns corresponds to the number of individual images and the number of rows corresponds to the number of y-coordinates for each function.  
25. **YX_Random_CDF_x** – subset of global x-coordinates (YX_Global_x) corresponding to random Y-->X CDF y-coordinates (YX_Random_CDF_y_Matrix). This is a single column vector used for all random YàX CDFs.   
26. **YX_Random_CDF_y_Matrix** – y-coordinates for all individual random Y-->X CDFs that are all defined over the same global x-coordinate scheme (YX_Random_CDF_x). Functions for individual images are stored in every column and y-coordinates are stored along rows. Thus, the number of columns corresponds to the number of individual images and the number of rows corresponds to the number of y-coordinates for each function.  
27. **YX_Delta_CDF_x** – subset of global x-coordinates (YX_Global_x) corresponding to Y-->X delta function y-coordinates (YX_Delta_CDF_y_Matrix). This is a single column vector used for all Y-->X delta functions.   
28. **YX_Delta_CDF_y_Matrix** – y-coordinates for all individual Y-->X delta functions that are all defined over the same global x-coordinate scheme (YX_Delta_CDF_x). Functions for individual images are stored in every column and y-coordinates are stored along rows. Thus, the number of columns corresponds to the number of individual images and the number of rows corresponds to the number of y-coordinates for each function.  
29. **YX_Observed_CDF_y_Median** – y-coordinates for the median Y-->X observed CDF. Computed as the weighted median for every row of YX_Observed_CDF_y_Matrix. Corresponding x-coordinates are stored in YX_Observed_CDF_x.  
30. **YX_Observed_CDF_y_Lower** – y-coordinates for the lower quantile of the median Y-->X observed CDF. Computed as the weighted lower quantile for every row of YX_Observed_CDF_y_Matrix. Corresponding x-coordinates are stored in YX_Observed_CDF_x.  
31. **YX_Observed_CDF_y_Upper** – y-coordinates for the upper quantile of the median Y-->X observed CDF. Computed as the weighted upper quantile for every row of YX_Observed_CDF_y_Matrix. Corresponding x-coordinates are stored in YX_Observed_CDF_x.  
32. **YX_Random_CDF_y_Median** – y-coordinates for the median Y-->X random CDF. Computed as the weighted median for every row of YX_Random_CDF_y_Matrix. Corresponding x-coordinates are stored in YX_Random_CDF_x.  
33. **YX_Random_CDF_y_Lower** – y-coordinates for the lower quantile of the median Y-->X random CDF. Computed as the weighted lower quantile for every row of YX_Random_CDF_y_Matrix. Corresponding x-coordinates are stored in YX_Random_CDF_x.  
34. **YX_Random_CDF_y_Upper** – y-coordinates for the upper quantile of the median Y-->X random CDF. Computed as the weighted upper quantile for every row of YX_Random_CDF_y_Matrix. Corresponding x-coordinates are stored in YX_Random_CDF_x.  
35. **YX_Delta_CDF_y_Median** – y-coordinates for the median Y-->X delta function. Computed as the weighted median for every row of YX_Delta_CDF_y_Matrix. Corresponding x-coordinates are stored in YX_Delta_CDF_x.  
36. **YX_Delta_CDF_y_Lower** – y-coordinates for the lower quantile of the median Y-->X delta function. Computed as the weighted lower quantile for every row of YX_Delta_CDF_y_Matrix. Corresponding x-coordinates are stored in YX_Delta_CDF_x.  
37. **YX_Delta_CDF_y_Upper** – y-coordinates for the upper quantile of the median Y-->X delta function. Computed as the weighted upper quantile for every row of YX_Delta_CDF_y_Matrix. Corresponding x-coordinates are stored in YX_Delta_CDF_x.  
38. **YX_Sample_Event_Count** – List of X-event counts for each individual image, used to generate all Y-->X functions. This is a column vector and each element corresponds to the columns in all Y-->X function matrices.   
39. **YX_Sample_Spatial_Association_Index** – List of Y-->X spatial association indices for individual images. Computed as the absolute maximum of each column in YX_Delta_CDF_y_Matrix.   
40. **YX_Global_Spatial_Association_Index** – Global spatial association index of Y-->X given all images in the batch. Computed as the absolute maximum of YX_Delta_CDF_y_Median.  
41. **YX_Global_CSR_Verdict** – verdict of the test for CSR in Y-->X. Equal to ‘true’ if the Y-->X median delta function quantile envelope (YX_Delta_CDF_y_Lower and YX_Delta_CDF_y_Upper) ever deviates from 0, and equal to ‘false’ if the envelope overlaps with 0 for the full distance range.  
