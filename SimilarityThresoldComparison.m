clc; clearvars; close all

%% Load
load SimilarityNetwork_FDR.mat
nC    = numel(countries);
alpha = 0.01;

%% ============================================================
%% Ensure diagonals are zero
S(1:nC+1:end) = 0;

%% ============================================================
%% Rebuild all filtered networks from S and P consistently
%% Original: threshold at P < alpha (uncorrected)
S_orig             = S;
S_orig(P >= alpha) = 0;
S_orig(1:nC+1:end) = 0;

%% FDR (Benjamini-Hochberg)
ut   = triu(true(nC), 1);
idx  = find(ut);
pvec = P(idx);

[p_sorted, ~] = sort(pvec);
m = numel(p_sorted);
bh_thresh = (1:m)'/m * alpha;
k = find(p_sorted <= bh_thresh, 1, 'last');

if isempty(k)
    p_cutoff = 0;
else
    p_cutoff = p_sorted(k);
end

S_fdr              = S;
S_fdr(P > p_cutoff) = 0;
S_fdr(1:nC+1:end)  = 0;

%% FWER (Bonferroni)
alpha_bonf = alpha / m;

S_fwer                  = S;
S_fwer(P >= alpha_bonf) = 0;
S_fwer(1:nC+1:end)     = 0;

%% ============================================================
%% Multiple Quantile-based pruning levels
%% Define quantile levels (proportion of edges to REMOVE)
quantile_levels = [0.01, 0.05, 0.10, 0.25, 0.50];
nQuantiles = numel(quantile_levels);

D_ut = D(idx);

%% Create cell arrays for quantile-pruned networks
S_quantile = cell(nQuantiles, 1);
A_quantile = cell(nQuantiles, 1);
quant_thresholds = zeros(nQuantiles, 1);

for qq = 1:nQuantiles
    quant_thresholds(qq) = quantile(D_ut, quantile_levels(qq));
    
    S_temp = S;
    S_temp(D >= quant_thresholds(qq)) = 0;
    S_temp(1:nC+1:end) = 0;
    S_quantile{qq} = S_temp;
    A_quantile{qq} = S_temp > 0;
end

%% ============================================================
%% Binary adjacency matrices for statistical methods
A_orig     = S_orig > 0;
A_fdr      = S_fdr > 0;
A_fwer     = S_fwer > 0;

%% ============================================================
%% Combine all methods
%% Statistical methods first, then quantile methods
stat_methods = {'Original', 'FDR', 'FWER'};
quant_methods = cell(nQuantiles, 1);
for qq = 1:nQuantiles
    quant_methods{qq} = sprintf('Quant_%d%%', round(quantile_levels(qq)*100));
end

methods  = [stat_methods, quant_methods'];
S_all    = [{S_orig, S_fdr, S_fwer}, S_quantile'];
A_all    = [{A_orig, A_fdr, A_fwer}, A_quantile'];
nMethods = numel(methods);

%% ============================================================
%% Compute metrics for all methods
density      = zeros(nMethods, 1);
total_weight = zeros(nMethods, 1);
avg_weight   = zeros(nMethods, 1);
strength     = zeros(nC, nMethods);
r_str        = zeros(nMethods, 1);
p_str        = zeros(nMethods, 1);
jaccard      = zeros(nMethods, 1);
r_wt         = zeros(nMethods, 1);
p_wt         = zeros(nMethods, 1);
fro_dist     = zeros(nMethods, 1);

nPossible = nnz(ut);   % total possible edges

for kk = 1:nMethods
    Sk = S_all{kk};
    Ak = A_all{kk};

    %% Edge density
    density(kk) = nnz(Ak & ut) / nPossible;

    %% Total weight
    total_weight(kk) = sum(Sk(ut));

    %% Average weight
    mask = Ak & ut;
    if nnz(mask) > 0
        avg_weight(kk) = mean(Sk(mask));
    else
        avg_weight(kk) = 0;
    end

    %% Node strength
    strength(:, kk) = sum(Sk, 2);

    %% Strength correlation vs Original
    if kk == 1
        r_str(kk) = 1; p_str(kk) = 0;
    elseif any(strength(:, kk))
        [r_str(kk), p_str(kk)] = corr(strength(:,1), strength(:,kk), ...
                                       'type', 'Spearman');
    else
        r_str(kk) = NaN; p_str(kk) = NaN;
    end

    %% Edge overlap (Jaccard vs Original)
    if kk == 1
        jaccard(kk) = 1;
    else
        inter_kk = nnz(A_orig & Ak & ut);
        union_kk = nnz((A_orig | Ak) & ut);
        if union_kk > 0
            jaccard(kk) = inter_kk / union_kk;
        else
            jaccard(kk) = 0;
        end
    end

    %% Weight correlation on common edges vs Original
    if kk == 1
        r_wt(kk) = 1; p_wt(kk) = 0;
    else
        common_kk  = (A_orig & Ak & ut);
        w_orig_kk  = S_orig(common_kk);
        w_method_kk = Sk(common_kk);
        if numel(w_orig_kk) > 2
            [r_wt(kk), p_wt(kk)] = corr(w_orig_kk, w_method_kk, ...
                                         'type', 'Spearman');
        else
            r_wt(kk) = NaN; p_wt(kk) = NaN;
        end
    end

    %% Frobenius distance vs Original
    fro_dist(kk) = norm(Sk - S_orig, 'fro');
end

%% ============================================================
%% Display all results
%% ============================================================

fprintf('\n===== THRESHOLD VALUES =====\n');
fprintf('FDR cutoff:        %.4g\n', p_cutoff);
fprintf('Bonferroni cutoff: %.4g\n', alpha_bonf);
fprintf('\nQuantile thresholds (distance):\n');
for qq = 1:nQuantiles
    fprintf('  %s: D < %.4f\n', quant_methods{qq}, quant_thresholds(qq));
end

fprintf('\n=== Full network comparison ===\n');

fprintf('\n--- Edge counts ---\n');
for kk = 1:nMethods
    fprintf('%-12s %d / %d\n', [methods{kk} ':'], nnz(A_all{kk} & ut), nPossible);
end

fprintf('\n--- Edge density ---\n');
for kk = 1:nMethods
    fprintf('%-12s %.3f\n', [methods{kk} ':'], density(kk));
end

fprintf('\n--- Total weight ---\n');
for kk = 1:nMethods
    fprintf('%-12s %.3f\n', [methods{kk} ':'], total_weight(kk));
end

fprintf('\n--- Average weight ---\n');
for kk = 1:nMethods
    fprintf('%-12s %.3f\n', [methods{kk} ':'], avg_weight(kk));
end

fprintf('\n--- Strength correlation (vs Original, Spearman) ---\n');
for kk = 2:nMethods
    fprintf('%-12s r = %.3f, p = %.3g\n', [methods{kk} ':'], r_str(kk), p_str(kk));
end

fprintf('\n--- Edge overlap (Jaccard vs Original) ---\n');
for kk = 2:nMethods
    fprintf('%-12s %.3f\n', [methods{kk} ':'], jaccard(kk));
end

fprintf('\n--- Weight correlation on common edges (vs Original, Spearman) ---\n');
for kk = 2:nMethods
    fprintf('%-12s r = %.3f, p = %.3g\n', [methods{kk} ':'], r_wt(kk), p_wt(kk));
end

fprintf('\n--- Frobenius distance (vs Original) ---\n');
for kk = 2:nMethods
    fprintf('%-12s %.3f\n', [methods{kk} ':'], fro_dist(kk));
end

%% ============================================================
%% Summary Table
%% ============================================================
fprintf('\n\n===== SUMMARY TABLE =====\n');
fprintf('%-12s %8s %8s %10s %10s %8s %10s\n', ...
    'Method', 'Edges', 'Density', 'TotWeight', 'AvgWeight', 'Jaccard', 'FrobDist');
fprintf('%s\n', repmat('-', 1, 70));
for kk = 1:nMethods
    fprintf('%-12s %8d %8.3f %10.3f %10.3f %8.3f %10.3f\n', ...
        methods{kk}, nnz(A_all{kk} & ut), density(kk), ...
        total_weight(kk), avg_weight(kk), jaccard(kk), fro_dist(kk));
end

%% ============================================================
%% Visualization
%% ============================================================
figure('Position', [100, 100, 1200, 400]);

%% Plot 1: Edge density comparison
subplot(1, 3, 1);
bar(density);
set(gca, 'XTickLabel', methods, 'XTickLabelRotation', 45);
ylabel('Edge Density');
title('Network Density by Pruning Method');
ylim([0, 1]);
grid on;

%% Plot 2: Jaccard similarity to Original
subplot(1, 3, 2);
bar(jaccard);
set(gca, 'XTickLabel', methods, 'XTickLabelRotation', 45);
ylabel('Jaccard Index');
title('Edge Overlap with Original Network');
ylim([0, 1]);
grid on;

%% Plot 3: Average edge weight
subplot(1, 3, 3);
bar(avg_weight);
set(gca, 'XTickLabel', methods, 'XTickLabelRotation', 45);
ylabel('Average Edge Weight');
title('Average Similarity of Retained Edges');
grid on;

sgtitle('Comparison of Graph Pruning Methods', 'FontSize', 14, 'FontWeight', 'bold');