clc;clearvars;close all
load('SimilarityNetwork_FDR.mat','Errore');
%% --- 1. Extract and Aggregate All Error Curves ---

% Check if Errore exists and is a cell array
if ~exist('Errore', 'var') || ~iscell(Errore)
    error('The variable "Errore" does not exist or is not a cell array.');
end

[num_rows, num_cols] = size(Errore);
all_error_curves = []; % Initialize a matrix to store all error curves as columns

fprintf('Extracting error curves from the Errore matrix...\n');
for i = 1:num_rows
    for j = 1:num_cols
        % Check if the cell is not empty (it should contain a 400x2 matrix)
        if ~isempty(Errore{i, j})
            % Extract the second column (the error values)
            current_error_curve = (Errore{i, j}(:, 2)+Errore{i, j}(:, 1))/2;
            
            % Append it as a new column to our collection
            all_error_curves = [all_error_curves, current_error_curve];
        end
    end
end

% Check if any data was found
if isempty(all_error_curves)
    error('No error data found in the "Errore" cell matrix. Is it populated correctly?');
end

num_curves = size(all_error_curves, 2);
fprintf('%d error curves found and extracted.\n', num_curves);

%% --- 2. Calculate Average and Standard Deviation ---

% Calculate the mean across all curves (dimension 2) for each iteration (row)
mean_error = mean(all_error_curves, 2);

% Calculate the standard deviation
std_error = std(all_error_curves, 0, 2);

% Get the number of iterations from the data
num_iterations = size(mean_error, 1);
iterations = 1:num_iterations;

%% --- 3. Plot the Average Error with Shaded Standard Deviation ---

fprintf('Generating the plot...\n');
figure('Color', 'white');
hold on;

% --- Create the shaded area for the standard deviation ---
% Define the x-coordinates for the fill polygon (forward and then backward)
fill_x = [iterations, fliplr(iterations)];
% Define the y-coordinates for the fill polygon (upper bound and then lower bound)
fill_y = [mean_error' + std_error', fliplr(mean_error' - std_error')];

% Use 'fill' to create the shaded region
h = fill(fill_x, fill_y, [0.8 0.8 1.0]); % A light blue color
set(h, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); % Make it transparent with no border

% --- Plot the mean error line on top ---
plot(iterations, mean_error, 'b-', 'LineWidth', 2);

% --- Formatting the Plot ---
title('Average Error Decrease Over Iterations', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('Iteration Number', 'FontSize', 12);
ylabel('Average Error', 'FontSize', 12);
legend('Standard Deviation', 'Mean Error', 'Location', 'northeast');
grid on;
box on;
axis tight; % Fit axes to the data
xlim([1 num_iterations]); % Ensure x-axis starts at 1

hold off;