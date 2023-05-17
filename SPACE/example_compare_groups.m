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
% Last Updated: 05/11/2023

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

% Independent (ind)
load([example_data, 'X_mask_list_ind.mat']);
load([example_data, 'Y_mask_list_ind.mat']);
load([example_data, 'ROI_mask_list_ind.mat']);

% Aggregated (agg)
load([example_data, 'X_mask_list_agg.mat']);
load([example_data, 'Y_mask_list_agg.mat']);
load([example_data, 'ROI_mask_list_agg.mat']);

%% Batch image SPACE analysis for each group of example images

[Batch_Results_dis, Single_Results_dis] = SPACE(X_mask_list_dis, Y_mask_list_dis, ROI_mask_list_dis);
[Batch_Results_ind, Single_Results_ind] = SPACE(X_mask_list_ind, Y_mask_list_ind, ROI_mask_list_ind);
[Batch_Results_agg, Single_Results_agg] = SPACE(X_mask_list_agg, Y_mask_list_agg, ROI_mask_list_agg);

%% Compare spatial association indices accross experimental groups using a weighted, 2-sided Student's T-test (ttest2w)

% dispersed VS independent
XY_p_dis_ind = ttest2w(Batch_Results_dis.XY_Sample_Spatial_Association_Index{1}, Batch_Results_ind.XY_Sample_Spatial_Association_Index{1}, Batch_Results_dis.XY_Sample_Event_Count{1}, Batch_Results_ind.XY_Sample_Event_Count{1}); % X-->Y 
YX_p_dis_ind = ttest2w(Batch_Results_dis.YX_Sample_Spatial_Association_Index{1}, Batch_Results_ind.YX_Sample_Spatial_Association_Index{1}, Batch_Results_dis.YX_Sample_Event_Count{1}, Batch_Results_ind.YX_Sample_Event_Count{1}); % Y-->X

% dispersed VS aggregated
XY_p_dis_agg = ttest2w(Batch_Results_dis.XY_Sample_Spatial_Association_Index{1}, Batch_Results_agg.XY_Sample_Spatial_Association_Index{1}, Batch_Results_dis.XY_Sample_Event_Count{1}, Batch_Results_agg.XY_Sample_Event_Count{1}); % X-->Y 
YX_p_dis_agg = ttest2w(Batch_Results_dis.YX_Sample_Spatial_Association_Index{1}, Batch_Results_agg.YX_Sample_Spatial_Association_Index{1}, Batch_Results_dis.YX_Sample_Event_Count{1}, Batch_Results_agg.YX_Sample_Event_Count{1}); % Y-->X

% independent VS aggregated
XY_p_ind_agg = ttest2w(Batch_Results_ind.XY_Sample_Spatial_Association_Index{1}, Batch_Results_agg.XY_Sample_Spatial_Association_Index{1}, Batch_Results_ind.XY_Sample_Event_Count{1}, Batch_Results_agg.XY_Sample_Event_Count{1}); % X-->Y 
YX_p_ind_agg = ttest2w(Batch_Results_ind.YX_Sample_Spatial_Association_Index{1}, Batch_Results_agg.YX_Sample_Spatial_Association_Index{1}, Batch_Results_ind.YX_Sample_Event_Count{1}, Batch_Results_agg.YX_Sample_Event_Count{1}); % Y-->X

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
xmax = max([Batch_Results_dis.XY_Global_x{1}(end), Batch_Results_dis.YX_Global_x{1}(end),...
            Batch_Results_ind.XY_Global_x{1}(end), Batch_Results_ind.YX_Global_x{1}(end),...
            Batch_Results_agg.XY_Global_x{1}(end), Batch_Results_agg.YX_Global_x{1}(end)]);

% Dispersed CDFs
subplot(2, 3, 1);
hold on;
[l(1), e(1)] = plot_error(Batch_Results_dis.XY_Observed_CDF_x{1}, Batch_Results_dis.XY_Observed_CDF_y_Median{1}, Batch_Results_dis.XY_Observed_CDF_y_Lower{1}, Batch_Results_dis.XY_Observed_CDF_y_Upper{1}); % dispersed observed CDF
[l(2), e(2)] = plot_error(Batch_Results_dis.XY_Random_CDF_x{1}, Batch_Results_dis.XY_Random_CDF_y_Median{1}, Batch_Results_dis.XY_Random_CDF_y_Lower{1}, Batch_Results_dis.XY_Random_CDF_y_Upper{1}); % dispersed random CDF
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
[l(1), e(1)] = plot_error(Batch_Results_ind.XY_Observed_CDF_x{1}, Batch_Results_ind.XY_Observed_CDF_y_Median{1}, Batch_Results_ind.XY_Observed_CDF_y_Lower{1}, Batch_Results_ind.XY_Observed_CDF_y_Upper{1}); % independent observed CDF
[l(2), e(2)] = plot_error(Batch_Results_ind.XY_Random_CDF_x{1}, Batch_Results_ind.XY_Random_CDF_y_Median{1}, Batch_Results_ind.XY_Random_CDF_y_Lower{1}, Batch_Results_ind.XY_Random_CDF_y_Upper{1}); % independent random CDF
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
[l(1), e(1)] = plot_error(Batch_Results_agg.XY_Observed_CDF_x{1}, Batch_Results_agg.XY_Observed_CDF_y_Median{1}, Batch_Results_agg.XY_Observed_CDF_y_Lower{1}, Batch_Results_agg.XY_Observed_CDF_y_Upper{1}); % aggregated observed CDF
[l(2), e(2)] = plot_error(Batch_Results_agg.XY_Random_CDF_x{1}, Batch_Results_agg.XY_Random_CDF_y_Median{1}, Batch_Results_agg.XY_Random_CDF_y_Lower{1}, Batch_Results_agg.XY_Random_CDF_y_Upper{1}); % aggregated random CDF
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
[p, e] = plot_error(Batch_Results_dis.XY_Delta_CDF_x{1}, Batch_Results_dis.XY_Delta_CDF_y_Median{1}, Batch_Results_dis.XY_Delta_CDF_y_Lower{1}, Batch_Results_dis.XY_Delta_CDF_y_Upper{1}); % dispersed delta CDF
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
[p, e] = plot_error(Batch_Results_ind.XY_Delta_CDF_x{1}, Batch_Results_ind.XY_Delta_CDF_y_Median{1}, Batch_Results_ind.XY_Delta_CDF_y_Lower{1}, Batch_Results_ind.XY_Delta_CDF_y_Upper{1}); % indpendent delta CDF
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
[p, e] = plot_error(Batch_Results_agg.XY_Delta_CDF_x{1}, Batch_Results_agg.XY_Delta_CDF_y_Median{1}, Batch_Results_agg.XY_Delta_CDF_y_Lower{1}, Batch_Results_agg.XY_Delta_CDF_y_Upper{1}); % aggregated delta CDF
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
[l(1), e(1)] = plot_error(Batch_Results_dis.YX_Observed_CDF_x{1}, Batch_Results_dis.YX_Observed_CDF_y_Median{1}, Batch_Results_dis.YX_Observed_CDF_y_Lower{1}, Batch_Results_dis.YX_Observed_CDF_y_Upper{1}); % dispersed observed CDF
[l(2), e(2)] = plot_error(Batch_Results_dis.YX_Random_CDF_x{1}, Batch_Results_dis.YX_Random_CDF_y_Median{1}, Batch_Results_dis.YX_Random_CDF_y_Lower{1}, Batch_Results_dis.YX_Random_CDF_y_Upper{1}); % dispersed random CDF
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
[l(1), e(1)] = plot_error(Batch_Results_ind.YX_Observed_CDF_x{1}, Batch_Results_ind.YX_Observed_CDF_y_Median{1}, Batch_Results_ind.YX_Observed_CDF_y_Lower{1}, Batch_Results_ind.YX_Observed_CDF_y_Upper{1}); % independent observed CDF
[l(2), e(2)] = plot_error(Batch_Results_ind.YX_Random_CDF_x{1}, Batch_Results_ind.YX_Random_CDF_y_Median{1}, Batch_Results_ind.YX_Random_CDF_y_Lower{1}, Batch_Results_ind.YX_Random_CDF_y_Upper{1}); % independent random CDF
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
[l(1), e(1)] = plot_error(Batch_Results_agg.YX_Observed_CDF_x{1}, Batch_Results_agg.YX_Observed_CDF_y_Median{1}, Batch_Results_agg.YX_Observed_CDF_y_Lower{1}, Batch_Results_agg.YX_Observed_CDF_y_Upper{1}); % aggregated observed CDF
[l(2), e(2)] = plot_error(Batch_Results_agg.YX_Random_CDF_x{1}, Batch_Results_agg.YX_Random_CDF_y_Median{1}, Batch_Results_agg.YX_Random_CDF_y_Lower{1}, Batch_Results_agg.YX_Random_CDF_y_Upper{1}); % aggregated random CDF
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
[p, e] = plot_error(Batch_Results_dis.YX_Delta_CDF_x{1}, Batch_Results_dis.YX_Delta_CDF_y_Median{1}, Batch_Results_dis.YX_Delta_CDF_y_Lower{1}, Batch_Results_dis.YX_Delta_CDF_y_Upper{1}); % dispersed delta CDF
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
[p, e] = plot_error(Batch_Results_ind.YX_Delta_CDF_x{1}, Batch_Results_ind.YX_Delta_CDF_y_Median{1}, Batch_Results_ind.YX_Delta_CDF_y_Lower{1}, Batch_Results_ind.YX_Delta_CDF_y_Upper{1}); % indpendent delta CDF
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
[p, e] = plot_error(Batch_Results_agg.YX_Delta_CDF_x{1}, Batch_Results_agg.YX_Delta_CDF_y_Median{1}, Batch_Results_agg.YX_Delta_CDF_y_Lower{1}, Batch_Results_agg.YX_Delta_CDF_y_Upper{1}); % aggregated delta CDF
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
s(1) = scatter(Batch_Results_dis.XY_Sample_Spatial_Association_Index{1}, Batch_Results_dis.YX_Sample_Spatial_Association_Index{1},'filled','markeredgecolor','k'); % dispersed
s(2) = scatter(Batch_Results_ind.XY_Sample_Spatial_Association_Index{1}, Batch_Results_ind.YX_Sample_Spatial_Association_Index{1},'filled','markeredgecolor','k'); % independent
s(3) = scatter(Batch_Results_agg.XY_Sample_Spatial_Association_Index{1}, Batch_Results_agg.YX_Sample_Spatial_Association_Index{1},'filled','markeredgecolor','k'); % aggregated
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
