% Example use of SPACE to characterize and compare the spatial association
% of point patterns X and Y between experimental groups. Tn this example,
% there are 3 distinct groups, those with dispersed (dis) signals,
% independent (ind) signals, and aggregated (agg) signals. Patterns are
% batch analyzed with SPACE, their spatial association indices are compared
% across groups with a weighted 2-sided ttest (ttest2w). And finally,
% median distribution functions and the spatial association index phase
% diagram are plotted from both perspectives of the SPACE analysis (X-->Y
% and Y-->X).
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
% Copyright 2024, Andrew Michael Soltisz. All rights reserved.
% This source code is licensed under the BSD-3-Clause License found in the
% LICENSE.txt file in the root directory of this source tree.

%% Prepare environment

close all; 
clear all; 
clc;

%% Import example data

% Generate example data if it doesn't exist
example_data = 'Example Data';
if ~isfolder(example_data)
    generate_example_data;
end
example_data = [example_data,filesep];

% Dispersed (dis)
load([example_data, 'X_mask_list_dis.mat']);
load([example_data, 'Y_mask_list_dis.mat']);
load([example_data, 'ROI_mask_list_dis.mat']);
n_images_dis = numel(X_mask_list_dis);

% Independent (ind)
load([example_data, 'X_mask_list_ind.mat']);
load([example_data, 'Y_mask_list_ind.mat']);
load([example_data, 'ROI_mask_list_ind.mat']);
n_images_ind = numel(X_mask_list_ind);

% Aggregated (agg)
load([example_data, 'X_mask_list_agg.mat']);
load([example_data, 'Y_mask_list_agg.mat']);
load([example_data, 'ROI_mask_list_agg.mat']);
n_images_agg = numel(X_mask_list_agg);

%% Perform SPACE Analysis of Individual Image Pairs

% dispersed dataset
wb = waitbar(0, sprintf("Analyzing dispersed image: %i/%i", 0, n_images_dis));
for i_image = 1:n_images_dis
    [X_Results_dis(i_image,:), Y_Results_dis(i_image,:)] = SPACE(X_mask_list_dis{i_image}, Y_mask_list_dis{i_image}); % default ROI is full image
    waitbar(i_image/n_images_dis, wb, sprintf("Analyzing dispersed image: %i/%i", i_image, n_images_dis));
end
close(wb);

% independent dataset
wb = waitbar(0, sprintf("Analyzing independent image: %i/%i", 0, n_images_ind));
for i_image = 1:n_images_ind
    [X_Results_ind(i_image,:), Y_Results_ind(i_image,:)] = SPACE(X_mask_list_ind{i_image}, Y_mask_list_ind{i_image}); % default ROI is full image
    waitbar(i_image/n_images_ind, wb, sprintf("Analyzing independent image: %i/%i", i_image, n_images_ind));
end
close(wb);

% aggregated dataset
wb = waitbar(0, sprintf("Analyzing aggregated image: %i/%i", 0, n_images_agg));
for i_image = 1:n_images_agg
    [X_Results_agg(i_image,:), Y_Results_agg(i_image,:)] = SPACE(X_mask_list_agg{i_image}, Y_mask_list_agg{i_image}); % default ROI is full image
    waitbar(i_image/n_images_agg, wb, sprintf("Analyzing aggregated image: %i/%i", i_image, n_images_agg));
end
close(wb);

%% Peform Batch Analysis of SPACE Results

% dispersed dataset
X_Batch_Results_dis = SPACE_batch(X_Results_dis);
Y_Batch_Results_dis = SPACE_batch(Y_Results_dis);

% independent dataset
X_Batch_Results_ind = SPACE_batch(X_Results_ind);
Y_Batch_Results_ind = SPACE_batch(Y_Results_ind);

% aggregated dataset
X_Batch_Results_agg = SPACE_batch(X_Results_agg);
Y_Batch_Results_agg = SPACE_batch(Y_Results_agg);

%% Compare spatial association indices accross experimental groups using a weighted, 2-sided Student's T-test (ttest2w)

% dispersed VS independent
X_p_dis_ind = ttest2w(X_Batch_Results_dis.Sample_Spatial_Association_Index{1}, X_Batch_Results_ind.Sample_Spatial_Association_Index{1}, X_Batch_Results_dis.Sample_Event_Count{1}, X_Batch_Results_ind.Sample_Event_Count{1}); % X-->Y 
Y_p_dis_ind = ttest2w(Y_Batch_Results_dis.Sample_Spatial_Association_Index{1}, Y_Batch_Results_ind.Sample_Spatial_Association_Index{1}, Y_Batch_Results_dis.Sample_Event_Count{1}, Y_Batch_Results_ind.Sample_Event_Count{1}); % Y-->X

% dispersed VS aggregated
X_p_dis_agg = ttest2w(X_Batch_Results_dis.Sample_Spatial_Association_Index{1}, X_Batch_Results_agg.Sample_Spatial_Association_Index{1}, X_Batch_Results_dis.Sample_Event_Count{1}, X_Batch_Results_agg.Sample_Event_Count{1}); % X-->Y 
Y_p_dis_agg = ttest2w(Y_Batch_Results_dis.Sample_Spatial_Association_Index{1}, Y_Batch_Results_agg.Sample_Spatial_Association_Index{1}, Y_Batch_Results_dis.Sample_Event_Count{1}, Y_Batch_Results_agg.Sample_Event_Count{1}); % Y-->X

% independent VS aggregated
X_p_ind_agg = ttest2w(X_Batch_Results_ind.Sample_Spatial_Association_Index{1}, X_Batch_Results_agg.Sample_Spatial_Association_Index{1}, X_Batch_Results_ind.Sample_Event_Count{1}, X_Batch_Results_agg.Sample_Event_Count{1}); % X-->Y 
Y_p_ind_agg = ttest2w(Y_Batch_Results_ind.Sample_Spatial_Association_Index{1}, Y_Batch_Results_agg.Sample_Spatial_Association_Index{1}, Y_Batch_Results_ind.Sample_Event_Count{1}, Y_Batch_Results_agg.Sample_Event_Count{1}); % Y-->X

%% Plot X-->Y Functions

% show masks
figure('position',[680   582   865   296]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact')
nexttile;
mask1_color = gen_overlay(X_mask_list_dis{1},Y_mask_list_dis{1});
imshow(imresize(mask1_color,4,'nearest'));
title("Dispersed Signals");
nexttile;
mask2_color = gen_overlay(X_mask_list_ind{1},Y_mask_list_ind{1});
imshow(imresize(mask2_color,4,'nearest'));
title("Independent Signals");
nexttile;
mask3_color = gen_overlay(X_mask_list_agg{1},Y_mask_list_agg{1});
imshow(imresize(mask3_color,4,'nearest'));
title("Aggregated Signals");

figure;
X_color = 'r';
Y_color = 'g';
sgtitle("X\rightarrowY Spatial Relationship");
% ensure common x-limit for all plots
xmax = max([X_Batch_Results_dis.Global_x{1}(end), Y_Batch_Results_dis.Global_x{1}(end),...
            X_Batch_Results_ind.Global_x{1}(end), Y_Batch_Results_ind.Global_x{1}(end),...
            X_Batch_Results_agg.Global_x{1}(end), Y_Batch_Results_agg.Global_x{1}(end)]);

% Dispersed CDFs
subplot(2, 3, 1);
hold on;
[l(1), e(1)] = plot_error(X_Batch_Results_dis.Observed_CDF_x{1}, X_Batch_Results_dis.Observed_CDF_y_Median{1}, X_Batch_Results_dis.Observed_CDF_y_Lower{1}, X_Batch_Results_dis.Observed_CDF_y_Upper{1}); % dispersed observed CDF
[l(2), e(2)] = plot_error(X_Batch_Results_dis.Random_CDF_x{1}, X_Batch_Results_dis.Random_CDF_y_Median{1}, X_Batch_Results_dis.Random_CDF_y_Lower{1}, X_Batch_Results_dis.Random_CDF_y_Upper{1}); % dispersed random CDF
l(1).Color = X_color;
l(2).Color = 'k';
l(2).LineStyle = '--';
e(1).FaceColor = X_color;
e(2).FaceColor = 'k';
legend(l,["Observed","Random"],'location','southeast');
title("Dispersed");
xlabel("Distance from Y", 'FontAngle', 'italic');
xlim([0, xmax]);
ylabel("P", 'FontAngle', 'italic');
ylim([0, 1]);
yticks(0:0.2:1);
grid on;
hold off;

% Independent CDFs
subplot(2,3 , 2);
hold on;
[l(1), e(1)] = plot_error(X_Batch_Results_ind.Observed_CDF_x{1}, X_Batch_Results_ind.Observed_CDF_y_Median{1}, X_Batch_Results_ind.Observed_CDF_y_Lower{1}, X_Batch_Results_ind.Observed_CDF_y_Upper{1}); % independent observed CDF
[l(2), e(2)] = plot_error(X_Batch_Results_ind.Random_CDF_x{1}, X_Batch_Results_ind.Random_CDF_y_Median{1}, X_Batch_Results_ind.Random_CDF_y_Lower{1}, X_Batch_Results_ind.Random_CDF_y_Upper{1}); % independent random CDF
l(1).Color = X_color;
l(2).Color = 'k';
l(2).LineStyle = '--';
e(1).FaceColor = X_color;
e(2).FaceColor = 'k';
title("Independent");
xlabel("Distance from Y", 'FontAngle', 'italic');
xlim([0, xmax]);
ylabel("P", 'FontAngle', 'italic');
ylim([0, 1]);
yticks(0:0.2:1);
grid on;
hold off;

% Aggregated CDFs
subplot(2,3,3);
hold on;
[l(1), e(1)] = plot_error(X_Batch_Results_agg.Observed_CDF_x{1}, X_Batch_Results_agg.Observed_CDF_y_Median{1}, X_Batch_Results_agg.Observed_CDF_y_Lower{1}, X_Batch_Results_agg.Observed_CDF_y_Upper{1}); % aggregated observed CDF
[l(2), e(2)] = plot_error(X_Batch_Results_agg.Random_CDF_x{1}, X_Batch_Results_agg.Random_CDF_y_Median{1}, X_Batch_Results_agg.Random_CDF_y_Lower{1}, X_Batch_Results_agg.Random_CDF_y_Upper{1}); % aggregated random CDF
l(1).Color = X_color;
l(2).Color = 'k';
l(2).LineStyle = '--';
e(1).FaceColor = X_color;
e(2).FaceColor = 'k';
title("Aggregated"); 
xlabel("Distance from Y", 'FontAngle', 'italic');
xlim([0, xmax]);
ylabel("P", 'FontAngle', 'italic');
ylim([0, 1]);
yticks(0:0.2:1);
grid on;
hold off;

% Dispersed delta CDFs
subplot(2, 3, 4);
hold on;
plot([0, xmax], [0, 0], 'k'); % plot black line across x-axis to highlight y=0
[p, e] = plot_error(X_Batch_Results_dis.Delta_CDF_x{1}, X_Batch_Results_dis.Delta_CDF_y_Median{1}, X_Batch_Results_dis.Delta_CDF_y_Lower{1}, X_Batch_Results_dis.Delta_CDF_y_Upper{1}); % dispersed delta CDF
p.Color = X_color;
e.FaceColor = X_color;
xlabel("Distance from Y", 'FontAngle', 'italic');
xlim([0, xmax]);
ylabel("\DeltaP", 'FontAngle', 'italic');
ylim([-1, 1]);
yticks(-1:0.25:1);
grid on;
hold off;

% Independent delta CDFs
subplot(2, 3, 5);
hold on;
plot([0, xmax], [0, 0], 'k'); % plot black line across x-axis to highlight y=0
[p, e] = plot_error(X_Batch_Results_ind.Delta_CDF_x{1}, X_Batch_Results_ind.Delta_CDF_y_Median{1}, X_Batch_Results_ind.Delta_CDF_y_Lower{1}, X_Batch_Results_ind.Delta_CDF_y_Upper{1}); % indpendent delta CDF
p.Color = X_color;
e.FaceColor = X_color;
xlabel("Distance from Y", 'FontAngle', 'italic');
xlim([0, xmax]);
ylabel("\DeltaP", 'FontAngle', 'italic');
ylim([-1, 1]);
yticks(-1:0.25:1);
grid on;
hold off;

% Aggregated delta CDFs
subplot(2, 3, 6);
hold on;
plot([0, xmax], [0, 0], 'k'); % plot black line across x-axis to highlight y=0
[p, e] = plot_error(X_Batch_Results_agg.Delta_CDF_x{1}, X_Batch_Results_agg.Delta_CDF_y_Median{1}, X_Batch_Results_agg.Delta_CDF_y_Lower{1}, X_Batch_Results_agg.Delta_CDF_y_Upper{1}); % aggregated delta CDF
p.Color = X_color;
e.FaceColor = X_color;
xlabel("Distance from Y", 'FontAngle', 'italic');
xlim([0, xmax]);
ylabel("\DeltaP", 'FontAngle', 'italic');
ylim([-1, 1]);
yticks(-1:0.25:1);
grid on;
hold off;

%% Plot Y-->X Functions

figure;
sgtitle("Y\rightarrowX Spatial Relationship");

% Dispersed CDFs
subplot(2, 3, 1);
hold on;
[l(1), e(1)] = plot_error(Y_Batch_Results_dis.Observed_CDF_x{1}, Y_Batch_Results_dis.Observed_CDF_y_Median{1}, Y_Batch_Results_dis.Observed_CDF_y_Lower{1}, Y_Batch_Results_dis.Observed_CDF_y_Upper{1}); % dispersed observed CDF
[l(2), e(2)] = plot_error(Y_Batch_Results_dis.Random_CDF_x{1}, Y_Batch_Results_dis.Random_CDF_y_Median{1}, Y_Batch_Results_dis.Random_CDF_y_Lower{1}, Y_Batch_Results_dis.Random_CDF_y_Upper{1}); % dispersed random CDF
l(1).Color = Y_color;
l(2).Color = 'k';
l(2).LineStyle = '--';
e(1).FaceColor = Y_color;
e(2).FaceColor = 'k';
legend(l,["Observed","Random"],'location','southeast');
title("Dispersed");
xlabel("Distance from X", 'FontAngle', 'italic');
xlim([0, xmax]);
ylabel("P", 'FontAngle', 'italic');
ylim([0, 1]);
yticks(0:0.2:1);
grid on;
hold off;

% Independent CDFs
subplot(2,3 , 2);
hold on;
[l(1), e(1)] = plot_error(Y_Batch_Results_ind.Observed_CDF_x{1}, Y_Batch_Results_ind.Observed_CDF_y_Median{1}, Y_Batch_Results_ind.Observed_CDF_y_Lower{1}, Y_Batch_Results_ind.Observed_CDF_y_Upper{1}); % independent observed CDF
[l(2), e(2)] = plot_error(Y_Batch_Results_ind.Random_CDF_x{1}, Y_Batch_Results_ind.Random_CDF_y_Median{1}, Y_Batch_Results_ind.Random_CDF_y_Lower{1}, Y_Batch_Results_ind.Random_CDF_y_Upper{1}); % independent random CDF
l(1).Color = Y_color;
l(2).Color = 'k';
l(2).LineStyle = '--';
e(1).FaceColor = Y_color;
e(2).FaceColor = 'k';
title("Independent");
xlabel("Distance from X", 'FontAngle', 'italic');
xlim([0, xmax]);
ylabel("P", 'FontAngle', 'italic');
ylim([0, 1]);
yticks(0:0.2:1);
grid on;
hold off;

% Aggregated CDFs
subplot(2,3,3);
hold on;
[l(1), e(1)] = plot_error(Y_Batch_Results_agg.Observed_CDF_x{1}, Y_Batch_Results_agg.Observed_CDF_y_Median{1}, Y_Batch_Results_agg.Observed_CDF_y_Lower{1}, Y_Batch_Results_agg.Observed_CDF_y_Upper{1}); % aggregated observed CDF
[l(2), e(2)] = plot_error(Y_Batch_Results_agg.Random_CDF_x{1}, Y_Batch_Results_agg.Random_CDF_y_Median{1}, Y_Batch_Results_agg.Random_CDF_y_Lower{1}, Y_Batch_Results_agg.Random_CDF_y_Upper{1}); % aggregated random CDF
l(1).Color = Y_color;
l(2).Color = 'k';
l(2).LineStyle = '--';
e(1).FaceColor = Y_color;
e(2).FaceColor = 'k';
title("Aggregated"); 
xlabel("Distance from X", 'FontAngle', 'italic');
xlim([0, xmax]);
ylabel("P", 'FontAngle', 'italic');
ylim([0, 1]);
yticks(0:0.2:1);
grid on;
hold off;

% Dispersed delta CDFs
subplot(2, 3, 4);
hold on;
plot([0, xmax], [0, 0], 'k'); % plot black line across x-axis to highlight y=0
[p, e] = plot_error(Y_Batch_Results_dis.Delta_CDF_x{1}, Y_Batch_Results_dis.Delta_CDF_y_Median{1}, Y_Batch_Results_dis.Delta_CDF_y_Lower{1}, Y_Batch_Results_dis.Delta_CDF_y_Upper{1}); % dispersed delta CDF
p.Color = Y_color;
e.FaceColor = Y_color;
xlabel("Distance from X", 'FontAngle', 'italic');
xlim([0, xmax]);
ylabel("\DeltaP", 'FontAngle', 'italic');
ylim([-1, 1]);
yticks(-1:0.25:1);
grid on;
hold off;

% Independent delta CDFs
subplot(2, 3, 5);
hold on;
plot([0, xmax], [0, 0], 'k'); % plot black line across x-axis to highlight y=0
[p, e] = plot_error(Y_Batch_Results_ind.Delta_CDF_x{1}, Y_Batch_Results_ind.Delta_CDF_y_Median{1}, Y_Batch_Results_ind.Delta_CDF_y_Lower{1}, Y_Batch_Results_ind.Delta_CDF_y_Upper{1}); % indpendent delta CDF
p.Color = Y_color;
e.FaceColor = Y_color;
xlabel("Distance from X", 'FontAngle', 'italic');
xlim([0, xmax]);
ylabel("\DeltaP", 'FontAngle', 'italic');
ylim([-1, 1]);
yticks(-1:0.25:1);
grid on;
hold off;

% Aggregated delta CDFs
subplot(2, 3, 6);
hold on;
plot([0, xmax], [0, 0], 'k'); % plot black line across x-axis to highlight y=0
[p, e] = plot_error(Y_Batch_Results_agg.Delta_CDF_x{1}, Y_Batch_Results_agg.Delta_CDF_y_Median{1}, Y_Batch_Results_agg.Delta_CDF_y_Lower{1}, Y_Batch_Results_agg.Delta_CDF_y_Upper{1}); % aggregated delta CDF
p.Color = Y_color;
e.FaceColor = Y_color;
xlabel("Distance from X", 'FontAngle', 'italic');
xlim([0, xmax]);
ylabel("\DeltaP", 'FontAngle', 'italic');
ylim([-1, 1]);
yticks(-1:0.25:1);
grid on;
hold off;

%% Phase Diagram of Spatial Association Indices for each group

figure;
hold on;
plot([-1, 1], [0, 0], 'k'); % plot black line across x-axis to highlight y=0
plot([0, 0], [-1, 1], 'k'); % plot black line across x-axis to highlight x=0
s(1) = scatter(X_Batch_Results_dis.Sample_Spatial_Association_Index{1}, Y_Batch_Results_dis.Sample_Spatial_Association_Index{1},'filled','markeredgecolor','k'); % dispersed
s(2) = scatter(X_Batch_Results_ind.Sample_Spatial_Association_Index{1}, Y_Batch_Results_ind.Sample_Spatial_Association_Index{1},'filled','markeredgecolor','k'); % independent
s(3) = scatter(X_Batch_Results_agg.Sample_Spatial_Association_Index{1}, Y_Batch_Results_agg.Sample_Spatial_Association_Index{1},'filled','markeredgecolor','k'); % aggregated
legend(s,["Dispersed","Independent","Aggregated"],'location','best');
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
