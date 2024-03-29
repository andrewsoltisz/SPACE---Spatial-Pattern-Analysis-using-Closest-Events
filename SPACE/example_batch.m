% Example use of SPACE to characterize the spatial association between
% point patterns X and Y from a batch of images. Pattern events are
% identified by true pixels in their corresponding image masks. 
%
% To analyze your own set of images, simply pack the masks for each signal
% into their own cell arrays and input them into the SPACE function.
% Spatial association indices can be compared across groups using the
% ttest2w function and median functions can easily be plotted with the
% plot_error function.
%
% See readme documentation for a detailed description of the output table
% fields.
%
% Spatial Pattern Analysis using Closest Events (SPACE)
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

clear all;
close all;
clc;

%% Import example data

% Generate example data if it doesn't exist
example_data = 'Example Data';
if ~isfolder(example_data)
    generate_example_data;
end
example_data = [example_data,filesep];

% Each 2D image is stored in its own element of a cell array
load([example_data, 'X_mask_list.mat']);
load([example_data, 'Y_mask_list.mat']);
load([example_data, 'ROI_mask_list.mat']);
n_images = numel(X_mask_list);

%% Batch image SPACE analysis

% Different use cases are listed below, un-comment one at a time
[Single_Results, Batch_Results] = SPACE(X_mask_list, Y_mask_list); % default ROI is full image
% [Single_Results, Batch_Results] = SPACE(X_mask_list, Y_mask_list, ROI_mask_list); % specify ROI
% [Single_Results, Batch_Results] = SPACE(X_mask_list, Y_mask_list, [], repelem(0.2,n_images,1)); % no ROI but pixel size of 0.2um is specified, input ROI as empty matrix
% [Single_Results, Batch_Results] = SPACE(X_mask_list, Y_mask_list, ROI_mask_list, repelem(0.2,n_images,1)); % Both ROI and pixel size of 0.2um are specified

%% Plot Results

% show masks
f = figure('position',[680   582   865   296]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact')
nexttile;
mask1_color = gen_overlay(X_mask_list{1},Y_mask_list{1});
imshow(imresize(mask1_color,4,'nearest'));
nexttile;
mask2_color = gen_overlay(X_mask_list{2},Y_mask_list{2});
imshow(imresize(mask2_color,4,'nearest'));
nexttile;
mask3_color = gen_overlay(X_mask_list{3},Y_mask_list{3});
imshow(imresize(mask3_color,4,'nearest'));
sgtitle(f,"Representative Images with Dispersed Patterns");

figure;
sgtitle("Batch SPACE Results");
X_color = 'r';
Y_color = 'g';
xmax = max([Batch_Results.XY_Global_x{1}(end), Batch_Results.XY_Global_x{1}(end)]); % ensure common x-limit for all plots

% X-->Y CDFs
subplot(2, 2, 1);
hold on;
[l(1), e(1)] = plot_error(Batch_Results.XY_Observed_CDF_x{1}, Batch_Results.XY_Observed_CDF_y_Median{1}, Batch_Results.XY_Observed_CDF_y_Lower{1}, Batch_Results.XY_Observed_CDF_y_Upper{1});
[l(2), e(2)] = plot_error(Batch_Results.XY_Random_CDF_x{1}, Batch_Results.XY_Random_CDF_y_Median{1}, Batch_Results.XY_Random_CDF_y_Lower{1}, Batch_Results.XY_Random_CDF_y_Upper{1}); 
l(1).Color = X_color;
l(2).Color = 'k';
l(2).LineStyle = '--';
e(1).FaceColor = X_color;
e(2).FaceColor = 'k';
legend(l,["Observed","Random"],'location','southeast');
title("X\rightarrowY Median CDFs");
xlabel("Distance from Y", 'FontAngle', 'italic');
xlim([0, xmax]);
ylabel("P", 'FontAngle', 'italic');
ylim([0, 1]);
yticks(0:0.2:1);
grid on;
hold off;

% Y-->X CDFs
subplot(2, 2, 2);
hold on;
[l(1), e(1)] = plot_error(Batch_Results.YX_Observed_CDF_x{1}, Batch_Results.YX_Observed_CDF_y_Median{1}, Batch_Results.YX_Observed_CDF_y_Lower{1}, Batch_Results.YX_Observed_CDF_y_Upper{1}); 
[l(2), e(2)] = plot_error(Batch_Results.YX_Random_CDF_x{1}, Batch_Results.YX_Random_CDF_y_Median{1}, Batch_Results.YX_Random_CDF_y_Lower{1}, Batch_Results.YX_Random_CDF_y_Upper{1}); 
l(1).Color = Y_color;
l(2).Color = 'k';
l(2).LineStyle = '--';
e(1).FaceColor = Y_color;
e(2).FaceColor = 'k';
legend(l,["Observed","Random"],'location','southeast');
title("Y\rightarrowX Median CDFs");
xlabel("Distance from X", 'FontAngle', 'italic');
xlim([0, xmax]);
ylabel("P", 'FontAngle', 'italic');
ylim([0, 1]);
yticks(0:0.2:1);
grid on;
hold off;

% X-->Y Delta CDFs
subplot(2, 2, 3);
hold on;
plot([0, xmax], [0, 0], 'k'); % plot y=0 to highlight x-axis
[p, e] = plot_error(Batch_Results.XY_Delta_CDF_x{1}, Batch_Results.XY_Delta_CDF_y_Median{1}, Batch_Results.XY_Delta_CDF_y_Lower{1}, Batch_Results.XY_Delta_CDF_y_Upper{1});
p.Color = X_color;
e.FaceColor = X_color;
title("X\rightarrowY Median \DeltaCDF");
xlabel("Distance from Y", 'FontAngle', 'italic');
xlim([0, xmax]);
ylabel("\DeltaP", 'FontAngle', 'italic');
ylim([-1, 1]);
yticks(-1:0.25:1);
grid on;
hold off;

% Y-->X Delta CDFs
subplot(2, 2, 4);
hold on;
plot([0, xmax], [0, 0], 'k'); % plot y=0 to highlight x-axis
[p, e] = plot_error(Batch_Results.YX_Delta_CDF_x{1}, Batch_Results.YX_Delta_CDF_y_Median{1}, Batch_Results.YX_Delta_CDF_y_Lower{1}, Batch_Results.YX_Delta_CDF_y_Upper{1});
p.Color = Y_color;
e.FaceColor = Y_color;
title("Y\rightarrowX Median \DeltaCDF");
xlabel("Distance from X", 'FontAngle', 'italic');
xlim([0, xmax]);
ylabel("\DeltaP", 'FontAngle', 'italic');
ylim([-1, 1]);
yticks(-1:0.25:1);
grid on;
hold off;

%% Phase Diagram of Spatial Association Indices

figure;
hold on;
plot([-1, 1], [0, 0], 'k'); % plot black line across x-axis to highlight y=0
plot([0, 0], [-1, 1], 'k'); % plot black line across x-axis to highlight x=0
scatter(Batch_Results.XY_Sample_Spatial_Association_Index{1}, Batch_Results.YX_Sample_Spatial_Association_Index{1},'filled','markeredgecolor','k'); % dispersed
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
