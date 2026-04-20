clc; clearvars; close all

%% Load data
load ReadyData.mat   % contains 'pdfs' and 'countries'

%% Parameters
nC = numel(countries);
alpha_values = [0.01, 0.05, 0.1];       % Significance levels
epsilon_values = [1e-4, 1e-5, 1e-6];    % Sinkhorn regularization values
B = 500;           % number of permutations
tol = 1e-8;
niter = 1000;

%% Define 2D cost matrix (based on bin positions)
edges = 0:1:10;       
xgrid = edges(1:end-1) + 0.5;  
ygrid = edges(1:end-1) + 0.5;
[X, Y] = meshgrid(xgrid, ygrid);
XY = [X(:), Y(:)];
C = pdist2(XY, XY, 'euclidean');
C = C / max(C(:));   % scale to [0,1]

%% Initialize storage for results
results = struct();
result_counter = 1;

%% Main loops over parameters
for alpha_idx = 1:length(alpha_values)
    for epsilon_idx = 1:length(epsilon_values)
        
        qFDR = alpha_values(alpha_idx);
        epsilon = epsilon_values(epsilon_idx);
        
        fprintf('\n=== Running with alpha = %.3f, epsilon = %.1e ===\n', qFDR, epsilon);
        
        %% Initialize outputs for this parameter combination
        P = ones(nC);   % p-values
        D = zeros(nC);  % Wasserstein distances
        Errore = cell(nC);  % Similarities
        
        %% Main comparison loop
        tic;
        for i = 1:nC
            for j = i+1:nC
                p1 = pdfs{i}(:); 
                p2 = pdfs{j}(:);

                % Permutation test for similarity
                 [pval, W_obs,Erro] = permtest_sinkhorn_divergence( ...
                    p1, p2, C, epsilon, tol, niter, B);
        
                Errore{i,j}=Erro;
                P(i,j) = pval;
                P(j,i) = pval;
        
                D(i,j) = W_obs;
                D(j,i) = W_obs;
            end
        end
        elapsed_time = toc;
        fprintf('Pairwise comparisons completed in %.2f seconds\n', elapsed_time);
        
        %% Convert distances to similarities
        sigma = median(D(D>0));      % characteristic distance scale
        S = exp(-D.^2 / (2*sigma^2));
        
        %% ============================================================
        %% FDR (Benjamini–Hochberg)
        idx  = find(triu(ones(nC),1));
        pvec = P(idx);
        
        [p_sorted, order] = sort(pvec);
        m = numel(p_sorted);
        
        bh_thresh = (1:m)'/m * qFDR;
        k = find(p_sorted <= bh_thresh, 1, 'last');
        
        if isempty(k)
            p_cutoff = 0;
        else
            p_cutoff = p_sorted(k);
        end
        
        fprintf('FDR cutoff: %.4g\n', p_cutoff);

        %% Filter by statistical significance
        S_filtered = S;
        S_filtered(P >= p_cutoff) = 0;   % remove non-significant similarities
        
        
        %% Compute network properties
        % Create adjacency matrix (binary)
        A = S_filtered > 0;
        A(1:nC+1:end) = 0;  % remove self-loops
        
        % Number of nodes
        num_nodes = nC;
        
        % Number of links
        num_links = sum(A(:)) / 2;  % undirected, divide by 2
        
        % Network density
        max_possible_links = nC * (nC - 1) / 2;
        density = num_links / max_possible_links;
        
        % Total flow (sum of all similarity weights)
        total_flow = sum(S_filtered(:)) / 2;  % undirected
        
        % Average similarity (average weight of existing links)
        avg_similarity = mean(S_filtered(A > 0));
        
        % Store results
        results(result_counter).qFDR = qFDR;
        results(result_counter).epsilon = epsilon;
        results(result_counter).S_filtered = S_filtered;
        results(result_counter).P = P;
        results(result_counter).D = D;
        results(result_counter).Errore = Errore;
        
        % Network properties
        results(result_counter).num_nodes = num_nodes;
        results(result_counter).num_links = num_links;
        results(result_counter).density = density;
        results(result_counter).total_flow = total_flow;
        results(result_counter).avg_similarity = avg_similarity;
        
       
        % Average path length (using BFS since it's weighted)
        avg_path_length = calculate_avg_path_length(S_filtered);
        results(result_counter).avg_path_length = avg_path_length;
        
        % Clustering coefficient (weighted)
        clustering_coef = calculate_weighted_clustering(S_filtered);
        results(result_counter).clustering_coefficient = clustering_coef;
 
        
        %% Display results for this parameter combination
        fprintf('Network Properties:\n');
        fprintf('  Number of nodes: %d\n', results(result_counter).num_nodes);
        fprintf('  Number of links: %.0f\n', results(result_counter).num_links);
        fprintf('  Density: %.4f\n', results(result_counter).density);
        fprintf('  Total flow: %.4f\n', results(result_counter).total_flow);
        fprintf('  Average similarity: %.4f\n', results(result_counter).avg_similarity);
        fprintf('  Average path length: %.4f\n', results(result_counter).avg_path_length);
        fprintf('  Clustering coefficient: %.4f\n', results(result_counter).clustering_coefficient);
        
        %% Save individual results
        filename = sprintf('SimilarityNetwork_alpha%.3f_epsilon%.0e.mat', qFDR, epsilon);
        save(filename, 'S_filtered', 'P', 'D', 'countries', 'Errore', 'results');
        
        
        result_counter = result_counter + 1;
    end
end

%% Save all results
save('AllResults.mat', 'results', 'countries', 'alpha_values', 'epsilon_values');

%% Create summary table
summary_table = struct2table(results);
disp('Summary of all results:');
disp(summary_table(:, {'qFDR', 'epsilon', 'num_nodes', 'num_links', 'density', ...
    'total_flow', 'avg_similarity', 'avg_path_length', 'clustering_coefficient'}));


%% ========================================================================
%% Helper function: Calculate average path length
function avg_length = calculate_avg_path_length(W)
% Calculate average shortest path length for weighted network
% W: weighted adjacency matrix (similarities)

n = size(W, 1);
D = zeros(n);  % Distance matrix

% Convert similarities to distances (higher similarity = shorter distance)
% Add small epsilon to avoid division by zero
W_temp = W + eps;
D_temp = 1 ./ W_temp;
D_temp(W == 0) = Inf;  % No connection = infinite distance
D_temp(1:n+1:end) = 0;  % Distance to self = 0

% Floyd-Warshall algorithm for all-pairs shortest paths
D = D_temp;
for k = 1:n
    for i = 1:n
        for j = 1:n
            if D(i,k) < Inf && D(k,j) < Inf
                if D(i,j) > D(i,k) + D(k,j)
                    D(i,j) = D(i,k) + D(k,j);
                end
            end
        end
    end
end

% Calculate average path length (excluding self and disconnected nodes)
path_lengths = [];
for i = 1:n
    for j = i+1:n
        if D(i,j) < Inf
            path_lengths = [path_lengths; D(i,j)];
        end
    end
end

if ~isempty(path_lengths)
    avg_length = mean(path_lengths);
else
    avg_length = NaN;
end
end

%% ========================================================================
%% Helper function: Calculate weighted clustering coefficient
function C = calculate_weighted_clustering(W)
% Calculate weighted clustering coefficient
% Based on Onnela et al. (2005)

n = size(W, 1);
C_i = zeros(n, 1);

% Normalize weights by maximum weight in network
W_norm = W / max(W(:));

for i = 1:n
    % Find neighbors of node i
    neighbors = find(W_norm(i,:) > 0);
    k = length(neighbors);
    
    if k < 2
        C_i(i) = 0;
    else
        % Calculate geometric mean of triangles
        sum_triangles = 0;
        for idx1 = 1:k
            for idx2 = idx1+1:k
                j = neighbors(idx1);
                m = neighbors(idx2);
                if W_norm(j,m) > 0
                    % Triangle weight: geometric mean of edge weights
                    triangle_weight = (W_norm(i,j) * W_norm(i,m) * W_norm(j,m))^(1/3);
                    sum_triangles = sum_triangles + triangle_weight;
                end
            end
        end
        C_i(i) = (2 * sum_triangles) / (k * (k-1));
    end
end

C = mean(C_i);
end


%% FUNCTION: permutation test using Sinkhorn divergence
function [pval, W_obs,Erro] = permtest_sinkhorn_divergence( ...
    p1, p2, C, epsilon, tol, niter, B)

p1 = p1 / sum(p1);
p2 = p2 / sum(p2);

% ---- observed divergence ----
[W_obs,Erro] = sinkhorn_divergence(p1, p2, C, epsilon, tol, niter);

% ---- permutation ----
combined = [p1(:); p2(:)];
n1 = numel(p1);

W_perm = zeros(B,1);

parfor b = 1:B

    idx = randperm(numel(combined));

    perm1 = combined(idx(1:n1));
    perm2 = combined(idx(n1+1:end));

    perm1 = perm1 / sum(perm1);
    perm2 = perm2 / sum(perm2);

    W_perm(b) = sinkhorn_divergence( ...
        perm1, perm2, C, epsilon, tol, niter);

end

% small divergence => similarity
pval = (1 + sum(W_perm <= W_obs)) / (B + 1);

end

%% ============================================================
%% FUNCTION: Sinkhorn divergence
function [W,Err] = sinkhorn_divergence(p1, p2, C, epsilon, tol, niter)

% OT(p1,p2)
[~,~,~,Err,OT12] = Sinkhorn_OT(C, epsilon, p1, p2, tol, niter);

% OT(p1,p1)
[~,~,~,~,OT11] = Sinkhorn_OT(C, epsilon, p1, p1, tol, niter);

% OT(p2,p2)
[~,~,~,~,OT22] = Sinkhorn_OT(C, epsilon, p2, p2, tol, niter);

% Sinkhorn divergence
W = OT12 - 0.5*OT11 - 0.5*OT22;

end



%% Entropic OT
function [T,a,b,Err,disto] = Sinkhorn_OT(C,epsilon,p,q, tol,niter)
    %% gibbs kernel
    K = exp(-C/epsilon);
    K(K<1e-200)=1e-200; % Safe
    q(q==0)=1e-50; % Safe
    p(p==0)=1e-50; % Safe
    a = ones(size(p));
    Err=nan(niter,2);
    %% Main
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
