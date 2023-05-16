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
% Last Updated: 05/11/2023

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

X_mask = imread([example_data, 'X_mask.tif']);
Y_mask = imread([example_data, 'Y_mask.tif']);
ROI_mask = imread([example_data, 'ROI_mask.tif']);

%% Single Image SPACE Analysis

Results = SPACE(X_mask, Y_mask, ROI_mask);

%% Plot Results

% create 1 figure with subplots for each plot
figure;
sgtitle("Single Image SPACE Results");
X_color = 'r';
Y_color = 'g';
xmax = max([Results.XY_Delta_CDF_x{1}(end), Results.XY_Delta_CDF_x{1}(end)]); % ensure common x-limit for all plots

% X-->Y CDFs
subplot(2,2,1);
hold on
p(1) = plot(Results.XY_Observed_x{1}, Results.XY_Observed_CDF_y{1}, X_color);
p(2) = plot(Results.XY_Random_x{1}, Results.XY_Random_CDF_y{1}, X_color, 'linestyle', '--');
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
p(1) = plot(Results.YX_Observed_x{1}, Results.YX_Observed_CDF_y{1}, Y_color);
p(2) = plot(Results.YX_Random_x{1}, Results.YX_Random_CDF_y{1}, Y_color, 'linestyle', '--');
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
plot(Results.XY_Delta_CDF_x{1}, Results.XY_Delta_CDF_y{1}, X_color);
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
plot(Results.YX_Delta_CDF_x{1}, Results.YX_Delta_CDF_y{1}, Y_color);
title("Y\rightarrowX \DeltaCDF");
xlabel("Distance from X");
xlim([0, xmax]);
ylabel("\DeltaP");
ylim([-1, 1])
yticks(-1:0.25:1);
grid on;
hold off;
