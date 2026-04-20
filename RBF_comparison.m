clc; clearvars; close all

%% ============================================================
%% Load data
load ReadyData.mat   % contains 'pdfs' and 'countries'

%% ============================================================
%% Parameters
nC      = numel(countries);
epsilon = 1e-5;
tol     = 1e-5;
niter   = 200;

%% ============================================================
%% Cost matrix (2D grid)
edges = 0:1:10;
xgrid = edges(1:end-1) + 0.5;
ygrid = edges(1:end-1) + 0.5;

[X,Y] = meshgrid(xgrid,ygrid);
XY = [X(:), Y(:)];

C = pdist2(XY, XY, 'euclidean');
C = C / max(C(:));

%% ============================================================
%% Initialize storage
D_wass = zeros(nC);     % Wasserstein distances
D_eucl = zeros(nC);     % Euclidean distances

%% ============================================================
%% MAIN LOOP: Compute both distance types
fprintf('\n');
fprintf('================================================================\n');
fprintf('COMPUTING DISTANCES\n');
fprintf('================================================================\n\n');

for i = 1:nC
    fprintf('Processing country %d/%d: %s\n', i, nC, countries{i});
    
    for j = i+1:nC
        
        p1 = pdfs{i}(:);
        p2 = pdfs{j}(:);
        
        % Normalize
        p1 = p1 / sum(p1);
        p2 = p2 / sum(p2);
        
        % Wasserstein distance (Sinkhorn divergence)
        W_obs = sinkhorn_divergence(p1, p2, C, epsilon, tol, niter);
        D_wass(i,j) = W_obs;
        D_wass(j,i) = W_obs;
        
        % Euclidean distance between proportion vectors
        D_eucl(i,j) = norm(p1 - p2);
        D_eucl(j,i) = D_eucl(i,j);
        
    end
end

fprintf('\nDistance computation complete.\n');

%% ============================================================
%% CREATE THREE SIMILARITY MATRICES
%% ============================================================

fprintf('\n');
fprintf('================================================================\n');
fprintf('CREATING SIMILARITY MATRICES\n');
fprintf('================================================================\n\n');

% Method 1: RBF kernel on WASSERSTEIN distance
sigma_wass = median(D_wass(D_wass>0));
S_wass_RBF = exp( -D_wass.^2 ./ (2*sigma_wass^2) );

% Method 2: RBF kernel on EUCLIDEAN distance
sigma_eucl = median(D_eucl(D_eucl>0));
S_eucl_RBF = exp( -D_eucl.^2 ./ (2*sigma_eucl^2) );

% Method 3: Raw WASSERSTEIN distance as dissimilarity
S_wass_raw = 1 ./ (1 + D_wass);

% Zero out diagonals
S_wass_RBF(1:nC+1:end) = 0;
S_eucl_RBF(1:nC+1:end) = 0;
S_wass_raw(1:nC+1:end) = 0;

fprintf('Method 1: Wasserstein + RBF kernel (sigma = %.4f)\n', sigma_wass);
fprintf('Method 2: Euclidean + RBF kernel   (sigma = %.4f)\n', sigma_eucl);
fprintf('Method 3: Wasserstein raw          (1/(1+d) transform)\n');

%% ============================================================
%% COMPARISON ANALYSIS
%% ============================================================

ut = triu(true(nC),1);  % Upper triangle mask

fprintf('\n');
fprintf('================================================================\n');
fprintf('TABLE 1: DISTANCE MEASURE COMPARISON\n');
fprintf('================================================================\n\n');

d_wass_vec = D_wass(ut);
d_eucl_vec = D_eucl(ut);

% Distance statistics
dist_stats = table();
dist_stats.Metric = {'Wasserstein'; 'Euclidean'};
dist_stats.Mean = [mean(d_wass_vec); mean(d_eucl_vec)];
dist_stats.Median = [median(d_wass_vec); median(d_eucl_vec)];
dist_stats.Std = [std(d_wass_vec); std(d_eucl_vec)];
dist_stats.Min = [min(d_wass_vec); min(d_eucl_vec)];
dist_stats.Max = [max(d_wass_vec); max(d_eucl_vec)];
dist_stats.Range = [range(d_wass_vec); range(d_eucl_vec)];

disp(dist_stats);

% Distance correlation
[r_dist_pearson, p_dist_p] = corr(d_wass_vec, d_eucl_vec, 'type', 'Pearson');
[r_dist_spearman, p_dist_s] = corr(d_wass_vec, d_eucl_vec, 'type', 'Spearman');

fprintf('\nCorrelation between Wasserstein and Euclidean distances:\n');
corr_dist = table();
corr_dist.Method = {'Pearson'; 'Spearman'};
corr_dist.Correlation = [r_dist_pearson; r_dist_spearman];
corr_dist.P_value = [p_dist_p; p_dist_s];
disp(corr_dist);

%% ============================================================
fprintf('\n');
fprintf('================================================================\n');
fprintf('TABLE 2: SIMILARITY MEASURE COMPARISON\n');
fprintf('================================================================\n\n');

s_wass_rbf_vec = S_wass_RBF(ut);
s_eucl_rbf_vec = S_eucl_RBF(ut);
s_wass_raw_vec = S_wass_raw(ut);

% Similarity statistics
sim_stats = table();
sim_stats.Method = {'Wass+RBF'; 'Eucl+RBF'; 'Wass_Raw'};
sim_stats.Mean = [mean(s_wass_rbf_vec); mean(s_eucl_rbf_vec); mean(s_wass_raw_vec)];
sim_stats.Median = [median(s_wass_rbf_vec); median(s_eucl_rbf_vec); median(s_wass_raw_vec)];
sim_stats.Std = [std(s_wass_rbf_vec); std(s_eucl_rbf_vec); std(s_wass_raw_vec)];
sim_stats.Min = [min(s_wass_rbf_vec); min(s_eucl_rbf_vec); min(s_wass_raw_vec)];
sim_stats.Max = [max(s_wass_rbf_vec); max(s_eucl_rbf_vec); max(s_wass_raw_vec)];
sim_stats.Range = [range(s_wass_rbf_vec); range(s_eucl_rbf_vec); range(s_wass_raw_vec)];

disp(sim_stats);

%% ============================================================
fprintf('\n');
fprintf('================================================================\n');
fprintf('TABLE 3: PAIRWISE SIMILARITY CORRELATIONS\n');
fprintf('================================================================\n\n');

[r_we, p_we] = corr(s_wass_rbf_vec, s_eucl_rbf_vec, 'type', 'Spearman');
[r_wr, p_wr] = corr(s_wass_rbf_vec, s_wass_raw_vec, 'type', 'Spearman');
[r_er, p_er] = corr(s_eucl_rbf_vec, s_wass_raw_vec, 'type', 'Spearman');

corr_sim = table();
corr_sim.Comparison = {'Wass+RBF vs Eucl+RBF'; 'Wass+RBF vs Wass_Raw'; 'Eucl+RBF vs Wass_Raw'};
corr_sim.Spearman_rho = [r_we; r_wr; r_er];
corr_sim.P_value = [p_we; p_wr; p_er];

disp(corr_sim);

%% ============================================================
fprintf('\n');
fprintf('================================================================\n');
fprintf('TABLE 4: NODE STRENGTH ANALYSIS\n');
fprintf('================================================================\n\n');

str_wass_rbf = sum(S_wass_RBF, 2);
str_eucl_rbf = sum(S_eucl_RBF, 2);
str_wass_raw = sum(S_wass_raw, 2);

% Node strength statistics
strength_stats = table();
strength_stats.Method = {'Wass+RBF'; 'Eucl+RBF'; 'Wass_Raw'};
strength_stats.Mean = [mean(str_wass_rbf); mean(str_eucl_rbf); mean(str_wass_raw)];
strength_stats.Median = [median(str_wass_rbf); median(str_eucl_rbf); median(str_wass_raw)];
strength_stats.Std = [std(str_wass_rbf); std(str_eucl_rbf); std(str_wass_raw)];
strength_stats.Min = [min(str_wass_rbf); min(str_eucl_rbf); min(str_wass_raw)];
strength_stats.Max = [max(str_wass_rbf); max(str_eucl_rbf); max(str_wass_raw)];

disp(strength_stats);

% Node strength correlations
[r_str_we, p_str_we] = corr(str_wass_rbf, str_eucl_rbf, 'type', 'Spearman');
[r_str_wr, p_str_wr] = corr(str_wass_rbf, str_wass_raw, 'type', 'Spearman');
[r_str_er, p_str_er] = corr(str_eucl_rbf, str_wass_raw, 'type', 'Spearman');

fprintf('\nNode Strength Correlations:\n');
corr_strength = table();
corr_strength.Comparison = {'Wass+RBF vs Eucl+RBF'; 'Wass+RBF vs Wass_Raw'; 'Eucl+RBF vs Wass_Raw'};
corr_strength.Spearman_rho = [r_str_we; r_str_wr; r_str_er];
corr_strength.P_value = [p_str_we; p_str_wr; p_str_er];

disp(corr_strength);

%% ============================================================
fprintf('\n');
fprintf('================================================================\n');
fprintf('TABLE 5: TOP 10 COUNTRIES BY NODE STRENGTH (Wass+RBF)\n');
fprintf('================================================================\n\n');

[sorted_str, idx_sorted] = sort(str_wass_rbf, 'descend');

% Create table with proper cell array handling
top_countries_data = cell(10, 5);
for k = 1:10
    top_countries_data{k, 1} = k;
    top_countries_data{k, 2} = countries{idx_sorted(k)};
    top_countries_data{k, 3} = sorted_str(k);
    top_countries_data{k, 4} = str_eucl_rbf(idx_sorted(k));
    top_countries_data{k, 5} = str_wass_raw(idx_sorted(k));
end

top_countries = cell2table(top_countries_data, ...
    'VariableNames', {'Rank', 'Country', 'Wass_RBF', 'Eucl_RBF', 'Wass_Raw'});

disp(top_countries);

%% ============================================================
fprintf('\n');
fprintf('================================================================\n');
fprintf('TABLE 6: MOST DISCORDANT COUNTRY PAIRS\n');
fprintf('(Standardized distance difference between Wasserstein & Euclidean)\n');
fprintf('================================================================\n\n');

% Standardize both distance measures
d_wass_std = (d_wass_vec - mean(d_wass_vec)) / std(d_wass_vec);
d_eucl_std = (d_eucl_vec - mean(d_eucl_vec)) / std(d_eucl_vec);

diff = abs(d_wass_std - d_eucl_std);
[sorted_diff, idx_sorted] = sort(diff, 'descend');

idx_pairs = find(ut);

% Create table with proper cell array handling
discordant_data = cell(10, 8);
for k = 1:10
    pair_idx = idx_sorted(k);
    linear_idx = idx_pairs(pair_idx);
    [i, j] = ind2sub([nC, nC], linear_idx);
    
    discordant_data{k, 1} = k;
    discordant_data{k, 2} = countries{i};
    discordant_data{k, 3} = countries{j};
    discordant_data{k, 4} = D_wass(i,j);
    discordant_data{k, 5} = D_eucl(i,j);
    discordant_data{k, 6} = sorted_diff(k);
    discordant_data{k, 7} = d_wass_std(pair_idx);
    discordant_data{k, 8} = d_eucl_std(pair_idx);
end

discordant_pairs = cell2table(discordant_data, ...
    'VariableNames', {'Rank', 'Country1', 'Country2', 'Wass_Dist', ...
                      'Eucl_Dist', 'Std_Diff', 'Wass_StdDist', 'Eucl_StdDist'});

disp(discordant_pairs);

%% ============================================================
fprintf('\n');
fprintf('================================================================\n');
fprintf('TABLE 7: MOST CONCORDANT COUNTRY PAIRS\n');
fprintf('(Smallest standardized distance difference)\n');
fprintf('================================================================\n\n');

[sorted_diff_asc, idx_sorted_asc] = sort(diff, 'ascend');

% Create table with proper cell array handling
concordant_data = cell(10, 6);
for k = 1:10
    pair_idx = idx_sorted_asc(k);
    linear_idx = idx_pairs(pair_idx);
    [i, j] = ind2sub([nC, nC], linear_idx);
    
    concordant_data{k, 1} = k;
    concordant_data{k, 2} = countries{i};
    concordant_data{k, 3} = countries{j};
    concordant_data{k, 4} = D_wass(i,j);
    concordant_data{k, 5} = D_eucl(i,j);
    concordant_data{k, 6} = sorted_diff_asc(k);
end

concordant_pairs = cell2table(concordant_data, ...
    'VariableNames', {'Rank', 'Country1', 'Country2', 'Wass_Dist', 'Eucl_Dist', 'Std_Diff'});

disp(concordant_pairs);

%% ============================================================
fprintf('\n');
fprintf('================================================================\n');
fprintf('TABLE 8: SIMILARITY RANGE COMPARISON\n');
fprintf('(Effect of RBF kernel vs raw transformation)\n');
fprintf('================================================================\n\n');

range_comparison = table();
range_comparison.Method = {'Wass+RBF'; 'Wass_Raw'; 'Ratio (RBF/Raw)'};
range_comparison.Min_Similarity = [min(s_wass_rbf_vec); min(s_wass_raw_vec); ...
    min(s_wass_rbf_vec)/min(s_wass_raw_vec)];
range_comparison.Max_Similarity = [max(s_wass_rbf_vec); max(s_wass_raw_vec); ...
    max(s_wass_rbf_vec)/max(s_wass_raw_vec)];
range_comparison.Range = [range(s_wass_rbf_vec); range(s_wass_raw_vec); ...
    range(s_wass_rbf_vec)/range(s_wass_raw_vec)];
range_comparison.Dynamic_Range = [max(s_wass_rbf_vec)/min(s_wass_rbf_vec); ...
    max(s_wass_raw_vec)/min(s_wass_raw_vec); NaN];

disp(range_comparison);

%% ============================================================
fprintf('\n');
fprintf('================================================================\n');
fprintf('TABLE 9: PERCENTILE ANALYSIS\n');
fprintf('================================================================\n\n');

percentiles = [10 25 50 75 90];

percentile_table = table();
percentile_table.Percentile = percentiles';

for p = 1:length(percentiles)
    percentile_table.Wass_Dist(p) = prctile(d_wass_vec, percentiles(p));
    percentile_table.Eucl_Dist(p) = prctile(d_eucl_vec, percentiles(p));
    percentile_table.Wass_RBF(p) = prctile(s_wass_rbf_vec, percentiles(p));
    percentile_table.Eucl_RBF(p) = prctile(s_eucl_rbf_vec, percentiles(p));
    percentile_table.Wass_Raw(p) = prctile(s_wass_raw_vec, percentiles(p));
end

disp(percentile_table);

%% ============================================================
fprintf('\n');
fprintf('================================================================\n');
fprintf('SUMMARY OF KEY FINDINGS\n');
fprintf('================================================================\n\n');

fprintf('1. DISTANCE METRICS:\n');
fprintf('   - Correlation (Spearman): %.3f\n', r_dist_spearman);
if r_dist_spearman > 0.9
    fprintf('   → Strong agreement: metrics rank pairs similarly\n');
elseif r_dist_spearman > 0.7
    fprintf('   → Moderate agreement: some divergence in rankings\n');
else
    fprintf('   → Weak agreement: metrics capture different structures\n');
end

fprintf('\n2. SIMILARITY TRANSFORMATIONS:\n');
fprintf('   - Wass+RBF vs Eucl+RBF: %.3f\n', r_we);
fprintf('   - Wass+RBF vs Wass_Raw: %.3f\n', r_wr);
if r_wr < 0.95
    fprintf('   → RBF kernel substantially modifies similarity structure\n');
else
    fprintf('   → RBF kernel preserves similarity structure\n');
end

fprintf('\n3. NODE STRENGTH AGREEMENT:\n');
fprintf('   - Wass+RBF vs Eucl+RBF: %.3f\n', r_str_we);
if r_str_we > 0.9
    fprintf('   → Country rankings highly consistent across methods\n');
else
    fprintf('   → Country rankings differ between methods\n');
end

fprintf('\n4. DYNAMIC RANGE (max/min ratio):\n');
fprintf('   - Wass+RBF: %.2f\n', max(s_wass_rbf_vec)/min(s_wass_rbf_vec));
fprintf('   - Wass_Raw: %.2f\n', max(s_wass_raw_vec)/min(s_wass_raw_vec));
if (max(s_wass_rbf_vec)/min(s_wass_rbf_vec)) > (max(s_wass_raw_vec)/min(s_wass_raw_vec))
    fprintf('   → RBF kernel provides greater dynamic range\n');
else
    fprintf('   → RBF kernel provides smaller dynamic range\n');
end

fprintf('\n');
fprintf('================================================================\n\n');

%% ============================================================
%% MINIMAL VISUALIZATION (just 2 key plots)
%% ============================================================

% Plot 1: Distance comparison
figure('Position', [100 100 800 350]);

subplot(1,2,1)
scatter(d_wass_vec, d_eucl_vec, 40, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('Wasserstein Distance', 'FontSize', 11);
ylabel('Euclidean Distance', 'FontSize', 11);
title(sprintf('Distance Correlation (ρ = %.3f)', r_dist_spearman), 'FontSize', 12);
grid on; box on; axis square;

subplot(1,2,2)
scatter(s_wass_rbf_vec, s_eucl_rbf_vec, 40, 'filled', 'MarkerFaceAlpha', 0.6);
xlabel('Wasserstein + RBF', 'FontSize', 11);
ylabel('Euclidean + RBF', 'FontSize', 11);
title(sprintf('Similarity Correlation (ρ = %.3f)', r_we), 'FontSize', 12);
grid on; box on; axis square;
hold on;
plot([0 1], [0 1], 'k--', 'LineWidth', 1.5);

sgtitle('Wasserstein vs Euclidean Comparison', 'FontSize', 13, 'FontWeight', 'bold');

% Plot 2: Effect of RBF kernelization
figure('Position', [150 150 800 350]);

subplot(1,2,1)
scatter(D_wass(ut), S_wass_raw(ut), 40, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('Wasserstein Distance', 'FontSize', 11);
ylabel('Raw Similarity: 1/(1+d)', 'FontSize', 11);
title('Raw Transformation', 'FontSize', 12);
grid on; box on;

subplot(1,2,2)
scatter(D_wass(ut), S_wass_RBF(ut), 40, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('Wasserstein Distance', 'FontSize', 11);
ylabel('RBF Similarity: exp(-d²/2σ²)', 'FontSize', 11);
title('RBF Kernel Transformation', 'FontSize', 12);
grid on; box on;

sgtitle('Why Use RBF Kernel?', 'FontSize', 13, 'FontWeight', 'bold');

%% ============================================================
%% Save
%% ============================================================
save('DistanceMetric_Comparison.mat', ...
    'S_wass_RBF', 'S_eucl_RBF', 'S_wass_raw', ...
    'D_wass', 'D_eucl', 'countries', ...
    'sim_stats', 'corr_sim', 'strength_stats', 'discordant_pairs');

fprintf('Results saved to DistanceMetric_Comparison.mat\n\n');

%% ============================================================
%% FUNCTIONS
%% ============================================================

function W = sinkhorn_divergence(p1, p2, C, epsilon, tol, niter)
    [~,~,~,~,OT12] = Sinkhorn_OT(C, epsilon, p1, p2, tol, niter);
    [~,~,~,~,OT11] = Sinkhorn_OT(C, epsilon, p1, p1, tol, niter);
    [~,~,~,~,OT22] = Sinkhorn_OT(C, epsilon, p2, p2, tol, niter);
    
    W = OT12 - 0.5*OT11 - 0.5*OT22;
end

function [T,a,b,Err,disto] = Sinkhorn_OT(C,epsilon,p,q, tol,niter)
    K = exp(-C/epsilon);
    K(K<1e-200)=1e-200;
    q(q==0)=1e-50;
    p(p==0)=1e-50;
    a = ones(size(p));
    Err=nan(niter,2);
    
    for i=1:niter
        b = q ./ (K'*a);
        if nargout>3 || tol>0
            Err(i,1) = norm(a.*(K*b)-p, 1);
        end
        a = p ./ (K*b);
        if nargout>3 || tol>0
            Err(i,2) = norm(b.*(K'*a)-q, 1);
        end
        if max(Err(i,:))<tol
            break;
        end
    end
    
    T = diag(a)*K*diag(b);
    aus=T.*log(T);
    aus(isinf(aus))=0;
    aus(isnan(aus))=0;
    
    disto = sum(sum(C.*T)) +sum(sum(epsilon*aus));
    disto(disto<0)=0;
end