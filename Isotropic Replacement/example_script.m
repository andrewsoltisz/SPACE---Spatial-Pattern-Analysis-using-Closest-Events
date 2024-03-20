% Example script illustrating how to use the isotropic_replacement.m
% function to correct spatial anisotropy in images
%
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% GitHub: https://github.com/andrewsoltisz/SPACE---Spatial-Pattern-Analysis-using-Closest-Events
% Publication: https://doi.org/10.1093/mam/ozae022
% Last Updated: 10/05/2023
%
% Copyright (C) 2024, Andrew Michael Soltisz. All rights reserved.
% This source code is licensed under the BSD-3-Clause License found in the
% LICENSE.txt file in the root directory of this source tree.

%% Prepare environment

close all;
clear all;
clc;

%% Import data

% define calibration
image_calibration_original = [0.024, 0.024, 0.15]; % microns per pixel (y,x,z) or (rows,cols,frames)

% imort images
mask_channel1_original = tiffreadVolume("example_image_channel1.tif");
mask_channel2_original = tiffreadVolume("example_image_channel2.tif");
im_sz_original = size(mask_channel1_original);
ROI_mask_original = true(im_sz_original); % define region of interest (ROI) as entire image

% show original masks
% combined masks to single RGB image
combined_mask = zeros([im_sz_original,3],'uint8');
combined_mask(:,:,:,1) = mask_channel1_original * 255; 
combined_mask(:,:,:,2) = mask_channel2_original * 255;
figure;
sliceViewer(combined_mask);
title("Original Anisotropic Masks");

%% Correct spatial anisotropy (multiple images at once)

% package images into a cell array
masks_original = {mask_channel1_original; mask_channel2_original};

% correct 
[masks_isotropic, ROI_mask_isotropic, image_calibration_new] = isotropic_replacement(masks_original, ROI_mask_original, image_calibration_original);
mask_channel1_isotropic = masks_isotropic{1};
mask_channel2_isotropic = masks_isotropic{2};
im_sz_new = size(mask_channel1_isotropic);

% show corrected masks
% combined masks to single RGB image
combined_mask = zeros([im_sz_new,3],'uint8');
combined_mask(:,:,:,1) = mask_channel1_isotropic * 255; 
combined_mask(:,:,:,2) = mask_channel2_isotropic * 255;
figure;
sliceViewer(combined_mask);
title("Corrected Isotropic Masks");

%% Perform SPACE with corrected images

% !WARNING!: ROI mask must be used with SPACE if anisotropy is corrected to
% specify which pixels in the corrected mask were part of the original
% mask. Otherwise, additional pixels will be included in the analysis and
% your sample size will be artificially increased.

pixel_size_isotropic = image_calibration_new(1);
Results = SPACE(mask_channel1_isotropic, mask_channel2_isotropic, ROI_mask_isotropic, pixel_size_isotropic);

%% Plot Results

% create 1 figure with subplots for each plot
figure;
sgtitle("Single Image SPACE Results");
channel1_color = 'r';
channel2_color = 'g';
xmax = max([Results.XY_Delta_CDF_x{1}(end), Results.XY_Delta_CDF_x{1}(end)]); % ensure common x-limit for all plots

% X-->Y CDFs
subplot(2,2,1);
hold on
p(1) = plot(Results.XY_Observed_x{1}, Results.XY_Observed_CDF_y{1}, channel1_color);
p(2) = plot(Results.XY_Random_x{1}, Results.XY_Random_CDF_y{1}, channel1_color, 'linestyle', '--');
legend(p,["Observed","Random"],'location','southeast');
title("X\rightarrowY CDFs");
xlabel("Distance from Y");
xlim([0, xmax]);
ylabel("P");
ylim([0, 1])
yticks(0:0.2:1);
grid on;
hold off;

% X-->Y CDFs
subplot(2,2,2);
hold on
p(1) = plot(Results.YX_Observed_x{1}, Results.YX_Observed_CDF_y{1}, channel2_color);
p(2) = plot(Results.YX_Random_x{1}, Results.YX_Random_CDF_y{1}, channel2_color, 'linestyle', '--');
legend(p,["Observed","Random"],'location','southeast');
title("Y\rightarrowX CDFs");
xlabel("Distance from X");
xlim([0, xmax]);
ylabel("P");
ylim([0, 1])
yticks(0:0.2:1);
grid on;
hold off;

% X-->Y CDFs
subplot(2,2,3);
hold on
plot([0,xmax], [0,0], 'k'); % plot y=0 to highlight x-axis
plot(Results.XY_Delta_CDF_x{1}, Results.XY_Delta_CDF_y{1}, channel1_color);
title("X\rightarrowY \DeltaCDF");
xlabel("Distance from Y");
xlim([0, xmax]);
ylabel("\DeltaP");
ylim([-1, 1])
yticks(-1:0.25:1);
grid on;
hold off;

% X-->Y CDFs
subplot(2,2,4);
hold on
plot([0,xmax], [0,0], 'k'); % plot y=0 to highlight x-axis
plot(Results.YX_Delta_CDF_x{1}, Results.YX_Delta_CDF_y{1}, channel2_color);
title("Y\rightarrowX \DeltaCDF");
xlabel("Distance from X");
xlim([0, xmax]);
ylabel("\DeltaP");
ylim([-1, 1])
yticks(-1:0.25:1);
grid on;
hold off;
