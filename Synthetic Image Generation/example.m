% Example script demonstrating how to use provided functions to generate
% synthetic images with tuned spatial association between signals.
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% GitHub: https://github.com/andrewsoltisz/SPACE---Spatial-Pattern-Analysis-using-Closest-Events
% Publication: https://doi.org/10.1101/2023.05.17.541131
% Last Updated: 10/05/2023
%
% Copyright (C) 2023, Andrew Michael Soltisz. All rights reserved.
% This source code is licensed under the BSD-3-Clause License found in the
% LICENSE.txt file in the root directory of this source tree.

%% Prepare environment

close all; 
clear all; 
clc;

%% Define image parameters (change as needed)

im_sz = [1000, 1000]; % pixels
X_conc = 0.1; % percent of total image, dynamic range [0,1]
Y_conc = 0.1; % percent of total image, dynamic range [0,1]
S = 0.8; % Y-->X Spatial proximity, dynamic range [-1,1]

%% Generate synthetic image

[Y_im1, X_im] = gen_synthetic_masks(im_sz, X_conc, Y_conc, S);
overlay_im1 = gen_overlay(X_im, Y_im1);
imshow(overlay_im1); 

%% Generate new Y-image using old X-image

Y_im2 = gen_synthetic_masks(im_sz, X_conc, Y_conc, S, X_im);
overlay_im2 = gen_overlay(X_im, Y_im2);
imshow(overlay_im2); 

%% View results

figure;
subplot(1,2,1);
imshow(overlay_im1); 
subplot(1,2,2);
imshow(overlay_im2); 

