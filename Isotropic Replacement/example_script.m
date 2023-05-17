% Example script illustrating how to use the isotropic_replacement.m
% function to correct spatial anisotropy in images
%
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% Last Updated: 05/16/2023

%% Prepare environment

close all;
clear all;
clc;

%% Import data

% define calibration
image_calibration_original = [0.024, 0.024, 0.15]; % microns per pixel (y,x,z) or (rows,cols,frames)

% imort images
image_channel1_original = tiffreadVolume("example_image_channel1.tif");
image_channel2_original = tiffreadVolume("example_image_channel2.tif");
image_channel3_original = tiffreadVolume("example_image_channel3.tif");

%% Correct spatial anisotropy (single image)

% correct 
[image_channel1_isotropic, image_calibration_new] = isotropic_replacement(image_channel1_original, image_calibration_original);

%% Correct spatial anisotropy (multiple images at once)

% package images into a cell array
images_original = {image_channel1_original; image_channel2_original; image_channel3_original};

% correct 
[images_isotropic, image_calibration_new] = isotropic_replacement(images_original, image_calibration_original);
