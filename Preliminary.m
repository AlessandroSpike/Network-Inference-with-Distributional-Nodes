% --- Main Script to Generate Plots ---
clc;clearvars;close all
% --- 1. Load and Prepare Data ---
% Load the data file which should contain a cell array 'pdfs' and a cell
% array of strings 'countries'.
load ReadyData.mat

% Dynamically get the number of countries and levels from the loaded data
% This is more robust than hardcoding the numbers.
num_countries = length(countries);
if ~isempty(pdfs)
    num_immig_levels = size(pdfs{1}, 1);
else
    error('The cell array "pdfs" is empty.');
end

% --- 2. Process Data ---
% We sum over the 'economy satisfaction' rows to get the marginal
% probability for each of the immigration answers.

% Matrix to hold the processed data: one row per country, one col per answer
processed_data = zeros(num_countries, num_immig_levels);
for i = 1:num_countries
    % Sum along the first dimension (rows) to integrate out economy satisfaction
    processed_data(i, :) = sum(pdfs{i}, 1);
end

% --- 3. Sort Data ---
% The chart is sorted by the total proportion of positive responses.
% Based on the image's 4 blue and 6 red/orange bars, and your 10 levels,
% we'll make the following split:
% - Cols 1-4: Negative sentiment (4 levels)
% - Cols 5-10: Positive sentiment (6 levels)

% **FIXED**: Calculate positive sentiment using only the positive columns
total_positive_sentiment = sum(processed_data(:, 5:10), 2);

[~, sort_idx] = sort(total_positive_sentiment, 'descend');

% Apply the sorting to our data and country names
sorted_data = processed_data(sort_idx, :);
sorted_countries = countries(sort_idx);


% --- 4. PLOT 1: Diverging Stacked Bar Chart (like the image) ---

figure('Position', [100, 100, 900, 800]); % Increased height for 28 countries
ax1 = axes;
hold(ax1, 'on'); % Hold on to plot multiple bar series

% Define the color maps, sampled from the image
% Colors for Negative sentiment (from most negative to least)
colors_neg = [
    30, 28, 70;    % Darkest Navy
    47, 52, 114;   % Dark Blue
    80, 99, 171;   % Medium Blue
    134, 153, 209; % Light Blue
] / 255; % Normalize to 0-1 range

% Colors for Positive sentiment (from least positive to most)
colors_pos = [
    253, 221, 143; % Pale Yellow
    252, 192, 106; % Light Yellow-Orange
    246, 160, 76;  % Orange
    220, 113, 58;  % Burnt Orange
    178, 59, 46;   % Red
    125, 23, 34;   % Dark Red
] / 255;

% Separate data into negative and positive parts for plotting
% **FIXED**: The indices now correctly match the 10-level data split.
% Negative data (cols 1-4). We reverse the order for plotting from the center out.
data_neg = sorted_data(:, 4:-1:1); 
% Positive data (cols 5-10)
data_pos = sorted_data(:, 5:10);

% Plot the positive bars (go from 0 to the right)
h_pos = barh(data_pos, 'stacked');

% Plot the negative bars (go from 0 to the left by making them negative)
h_neg = barh(-data_neg, 'stacked');

% Apply the custom colors
for i = 1:size(colors_pos, 1)
    h_pos(i).FaceColor = colors_pos(i, :);
    h_pos(i).EdgeColor = 'none'; % Remove black edges for a cleaner look
end
for i = 1:size(colors_neg, 1)
    h_neg(i).FaceColor = colors_neg(i, :);
    h_neg(i).EdgeColor = 'none';
end

% --- Formatting the Plot ---
set(ax1, 'YDir', 'reverse'); % Put the first country at the top
set(ax1, 'YTick', 1:num_countries);
set(ax1, 'YTickLabel', sorted_countries);
set(ax1, 'YAxisLocation', 'left');
ax1.YAxis.FontSize = 11;
ax1.TickLength = [0 0]; % Remove tick marks

% Format X-axis to show percentages, symmetric around 0
max_val = max(abs(ax1.XLim));
limit = ceil(max_val * 10) / 10; % Round up to nearest 0.1
if limit == 0; limit = 0.1; end % Avoid limit being zero if data is all zero
ax1.XLim = [-limit, limit];
% Generate ticks dynamically based on the limit
tick_step = 0.1;
if limit > 0.5; tick_step = 0.2; end
xticks_vals = -limit:tick_step:limit;
ax1.XTick = xticks_vals;

% Create percentage labels (e.g., 40, 20, 0, 20, 40)
xtick_labels = string(abs(xticks_vals * 100));
set(ax1, 'XTickLabel', xtick_labels);

% Add faint grid lines and clean up the look
ax1.XGrid = 'on';
ax1.YGrid = 'on';
ax1.GridColor = [0.8 0.8 0.8];
ax1.GridAlpha = 0.5;
box(ax1, 'off'); % Remove the box outline

title('Immigration Opinions by Country', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('Percentage of Respondents (%)', 'FontSize', 12);
axis tight
hold(ax1, 'off');




% --- 2. Process Data ---
% We sum over the 'economy satisfaction' rows to get the marginal
% probability for each of the immigration answers.

% Matrix to hold the processed data: one row per country, one col per answer
processed_data = zeros(num_countries, num_immig_levels);
for i = 1:num_countries
    % Sum along the first dimension (rows) to integrate out economy satisfaction
    processed_data(i, :) = sum(pdfs{i}, 2);
end

% --- 3. Sort Data ---
% The chart is sorted by the total proportion of positive responses.
% Based on the image's 4 blue and 6 red/orange bars, and your 10 levels,
% we'll make the following split:
% - Cols 1-4: Negative sentiment (4 levels)
% - Cols 5-10: Positive sentiment (6 levels)

% **FIXED**: Calculate positive sentiment using only the positive columns
total_positive_sentiment = sum(processed_data(:, 5:10), 2);

[~, sort_idx] = sort(total_positive_sentiment, 'descend');

% Apply the sorting to our data and country names
sorted_data = processed_data(sort_idx, :);
sorted_countries = countries(sort_idx);


% --- 4. PLOT 1: Diverging Stacked Bar Chart (like the image) ---

figure('Position', [100, 100, 900, 800]); % Increased height for 28 countries
ax1 = axes;
hold(ax1, 'on'); % Hold on to plot multiple bar series

% Define the color maps, sampled from the image
% Colors for Negative sentiment (from most negative to least)
colors_neg = [
    30, 28, 70;    % Darkest Navy
    47, 52, 114;   % Dark Blue
    80, 99, 171;   % Medium Blue
    134, 153, 209; % Light Blue
] / 255; % Normalize to 0-1 range

% Colors for Positive sentiment (from least positive to most)
colors_pos = [
    253, 221, 143; % Pale Yellow
    252, 192, 106; % Light Yellow-Orange
    246, 160, 76;  % Orange
    220, 113, 58;  % Burnt Orange
    178, 59, 46;   % Red
    125, 23, 34;   % Dark Red
] / 255;

% Separate data into negative and positive parts for plotting
% **FIXED**: The indices now correctly match the 10-level data split.
% Negative data (cols 1-4). We reverse the order for plotting from the center out.
data_neg = sorted_data(:, 4:-1:1); 
% Positive data (cols 5-10)
data_pos = sorted_data(:, 5:10);

% Plot the positive bars (go from 0 to the right)
h_pos = barh(data_pos, 'stacked');

% Plot the negative bars (go from 0 to the left by making them negative)
h_neg = barh(-data_neg, 'stacked');

% Apply the custom colors
for i = 1:size(colors_pos, 1)
    h_pos(i).FaceColor = colors_pos(i, :);
    h_pos(i).EdgeColor = 'none'; % Remove black edges for a cleaner look
end
for i = 1:size(colors_neg, 1)
    h_neg(i).FaceColor = colors_neg(i, :);
    h_neg(i).EdgeColor = 'none';
end

% --- Formatting the Plot ---
set(ax1, 'YDir', 'reverse'); % Put the first country at the top
set(ax1, 'YTick', 1:num_countries);
set(ax1, 'YTickLabel', sorted_countries);
set(ax1, 'YAxisLocation', 'left');
ax1.YAxis.FontSize = 11;
ax1.TickLength = [0 0]; % Remove tick marks

% Format X-axis to show percentages, symmetric around 0
max_val = max(abs(ax1.XLim));
limit = ceil(max_val * 10) / 10; % Round up to nearest 0.1
if limit == 0; limit = 0.1; end % Avoid limit being zero if data is all zero
ax1.XLim = [-limit, limit];
% Generate ticks dynamically based on the limit
tick_step = 0.1;
if limit > 0.5; tick_step = 0.2; end
xticks_vals = -limit:tick_step:limit;
ax1.XTick = xticks_vals;

% Create percentage labels (e.g., 40, 20, 0, 20, 40)
xtick_labels = string(abs(xticks_vals * 100));
set(ax1, 'XTickLabel', xtick_labels);

% Add faint grid lines and clean up the look
ax1.XGrid = 'on';
ax1.YGrid = 'on';
ax1.GridColor = [0.8 0.8 0.8];
ax1.GridAlpha = 0.5;
box(ax1, 'off'); % Remove the box outline

title('Economic Opinions by Country', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('Percentage of Respondents (%)', 'FontSize', 12);
axis tight
hold(ax1, 'off');

