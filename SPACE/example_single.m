% Example use of SPACE to characterize the spatial association between
% point patterns X and Y from a single image. Pattern events are identified
% by true pixels in their corresponding image masks.
%
% To analyze your own image, simply input the masks into the SPACE function
% and the results will be output in the form of a table.
%
% See readme documentation for a detailed description of the output table
% fields.
%
% Spatial Pattern Analysis using Closest Events (SPACE)
% Author: Andrew M. Soltisz
% Email: andysoltisz@gmail.com
% GitHub: https://github.com/andrewsoltisz/SPACE---Spatial-Pattern-Analysis-using-Closest-Events
% Publication: https://doi.org/10.1093/mam/ozae022
% Last Updated: 03/27/2025
%
% Copyright 2025, Andrew Michael Soltisz. All rights reserved.
% This source code is licensed under the BSD-3-Clause License found in the
% LICENSE.txt file in the root directory of this source tree.

%% Prepare environment

clear all;
% close all;
clc;

%% Import example data

% Generate example data if it doesn't exist
example_data = 'Example Data';
if ~isfolder(example_data)
    generate_example_data;
end
example_data = [example_data,filesep];

X_mask = imread([example_data, 'X_mask.tif']);
Y_mask = imread([example_data, 'Y_mask.tif']);
ROI_mask = imread([example_data, 'ROI_mask.tif']);
pixelSize = 0.2; % micrometers

%% Single Image SPACE Analysis 

% perform analysis from the perspective of both event species
Results_X = SPACE(X_mask, Y_mask, ROI_mask, pixelSize);  
Results_Y = SPACE(Y_mask, X_mask, ROI_mask, pixelSize); 

%% Plot Results

% show masks
figure('position',[680   582   865   296]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact')
nexttile;
X_mask_color = gen_overlay(X_mask,false(size(X_mask)));
imshow(imresize(X_mask_color,4,'nearest'));
title("X Mask");
nexttile;
Y_mask_color = gen_overlay(false(size(X_mask)),Y_mask);
imshow(imresize(Y_mask_color,4,'nearest'));
title("Y Mask");
nexttile;
XY_mask_color = gen_overlay(X_mask,Y_mask);
imshow(imresize(XY_mask_color,4,'nearest'));
title("Combined Mask");
sgtitle("Image with Aggregated Patterns");

% create 1 figure with subplots for each plot
figure;
sgtitle("Single Image SPACE Results");
X_color = 'r';
Y_color = 'g';
xmax = max([Results_X.Delta_CDF_x{1}(end), Results_Y.Delta_CDF_x{1}(end)]); % ensure common x-limit for all plots

% X-->Y CDFs
subplot(2,2,1);
hold on
p(1) = plot(Results_X.Observed_x{1}, Results_X.Observed_CDF_y{1}, X_color);
p(2) = plot(Results_X.Random_x{1}, Results_X.Random_CDF_y{1}, X_color, 'linestyle', '--');
legend(p,["Observed","Random"],'location','southeast');
title("X\rightarrowY CDFs");
xlabel("Distance from Y");
xlim([0, xmax]);
ylabel("P");
ylim([0, 1])
yticks(0:0.2:1);
grid on;
hold off;

% Y-->Y CDFs
subplot(2,2,2);
hold on
p(1) = plot(Results_Y.Observed_x{1}, Results_Y.Observed_CDF_y{1}, Y_color);
p(2) = plot(Results_Y.Random_x{1}, Results_Y.Random_CDF_y{1}, Y_color, 'linestyle', '--');
legend(p,["Observed","Random"],'location','southeast');
title("Y\rightarrowX CDFs");
xlabel("Distance from X");
xlim([0, xmax]);
ylabel("P");
ylim([0, 1])
yticks(0:0.2:1);
grid on;
hold off;

% X-->Y Delta CDFs
subplot(2,2,3);
hold on
plot([0,xmax], [0,0], 'k'); % plot y=0 to highlight x-axis
plot(Results_X.Delta_CDF_x{1}, Results_X.Delta_CDF_y{1}, X_color);
title("X\rightarrowY \DeltaCDF");
xlabel("Distance from Y");
xlim([0, xmax]);
ylabel("\DeltaP");
ylim([-1, 1])
yticks(-1:0.25:1);
grid on;
hold off;

% Y-->Y Delta CDFs
subplot(2,2,4);
hold on
plot([0,xmax], [0,0], 'k'); % plot y=0 to highlight x-axis
plot(Results_Y.Delta_CDF_x{1}, Results_Y.Delta_CDF_y{1}, Y_color);
title("Y\rightarrowX \DeltaCDF");
xlabel("Distance from X");
xlim([0, xmax]);
ylabel("\DeltaP");
ylim([-1, 1])
yticks(-1:0.25:1);
grid on;
hold off;

%% Phase Diagram of Spatial Association Indices

figure;
hold on;
plot([-1, 1], [0, 0], 'k'); % plot black line across x-axis to highlight y=0
plot([0, 0], [-1, 1], 'k'); % plot black line across x-axis to highlight x=0
scatter(Results_X.Spatial_Association_Index,Results_Y.Spatial_Association_Index,100,'filled','markeredgecolor','k'); % aggregated
title("Spatial Association Index Phase Diagram");
xlabel("X\rightarrowY Spatial Association");
xlim([-1, 1]);
xticks(-1:0.25:1);
ylabel("Y\rightarrowX Spatial Association");
ylim([-1, 1]);
yticks(-1:0.25:1);
grid on;
axis square;
hold off;
