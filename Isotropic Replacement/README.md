## Authorship:  
**Author**: Andrew M. Soltisz  
**Contact**: andysoltisz@gmail.com  
**LinkedIn**: https://www.linkedin.com/in/andrew-soltisz/  
**Publication**: This work is described in more detail at https://doi.org/10.1101/2023.05.17.541131  

Copyright (C) 2023, Andrew Michael Soltisz.  
All rights reserved.   
The source code within this directory is licensed under the BSD-3-Clause License found in the LICENSE.txt file.


## Description:  
Distance transformations require spatially isotropic images to produce accurate Euclidean distances. However, due to the elongation of the point-spread-function in the axial dimension, it is common for volumetric microscopy modalities to have poorer resolution in the z-dimension. As a result, sampling is often coarser in the axial dimension (z) than in the lateral dimensions (x, y), rendering the resulting digital images spatially anisotropic. While this could be corrected by simply resampling the image, the interpolation step will alter both signal morphology and the number of pixels that compose the signals, and thus reduce statistical power of subsequent analyses. Instead, image spatial anisotropy can be corrected, while avoiding these issues, using a technique we call 'isotropic signal replacement'. Following image segmentation, a new blank image template is created with the same overall dimensions as the input image but with isotropic sampling, i.e. the z-dimension is sampled at the same spatial frequency as the lateral dimensions. Segmented pixels from the input image are then inserted back into this template image at pixel locations that are closest to the signal’s original Euclidean coordinates – in other words, pixels in the new image with centroid positions closest to those of the original signal-positive pixels in the raw image. It should be noted that positional rounding during the pixel replacement step can lead to small deviations in these pixel positions between the input image and the newly created image, thus altering measured distance. However, if the z-dimension is subdivided by a factor that is less than or equal to the lowest image resolution divided by root three, then the distance error between any two pixels within the new image will be less than the optical resolution of the image. Thus, any resulting errors in measured distances would fall below the least significant digit.
