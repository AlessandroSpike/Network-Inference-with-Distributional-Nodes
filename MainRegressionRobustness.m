clc; clearvars; close all;

%% Load data
load ReadyData.mat   % contains 'pdfs' and 'countries'

%% Parameters
nC = numel(countries);
alpha_values = [0.01, 0.05, 0.1];       % Significance thresholds
epsilon_values = [1e-4, 1e-5, 1e-6];    % Sinkhorn regularization values
alphav = 1;         % step size for barycenter regression
tol      = 1e-8;     % Sinkhorn tolerance
tol2     = 1e-8;     % gradient tolerance
maxiter  = 50;      % max iterations for barycenter regression
dxmin = 1e-8;       % minimum perturbation
theta = Inf;        % unused
L = 50;            % gradient iterations for barycenter
maxK = 27;          % maximum number of neighbors
B = 500;            % number of Monte Carlo samples

%% Define 2D cost matrix (bin positions)
edges = 0:1:10;
xgrid = edges(1:end-1) + 0.5;
ygrid = edges(1:end-1) + 0.5;
[X, Y] = meshgrid(xgrid, ygrid);
XY = [X(:), Y(:)];
C = pdist2(XY, XY, 'euclidean');
C = C / max(C(:));

%% Reshape data
threeD_Array = cat(3, pdfs{:});
pdfs2 = reshape(permute(threeD_Array, [3, 1, 2]), nC, []);
pdfs2 = pdfs2';  % columns = countries

%% Initialize storage for results
results = struct();
result_counter = 1;

%% Main loops over parameters
for alpha_idx = 1:length(alpha_values)
    for epsilon_idx = 1:length(epsilon_values)
        
        alpha = alpha_values(alpha_idx);
        epsilon = epsilon_values(epsilon_idx);
        
        fprintf('\n=== Running with alpha = %.3f, epsilon = %.1e ===\n', alpha, epsilon);
        
        %% Load or compute distance matrix for this epsilon
        % Try to load precomputed distances or compute new ones
        dist_filename = sprintf('SimilarityNetwork_alpha0.010_epsilon%.0e.mat', epsilon);
        load(dist_filename, 'D');
        fprintf('Loaded distance matrix from %s\n', dist_filename);
        
        %% Initialize outputs for this parameter combination
        AdjMatrix = zeros(nC);
        Norma1 = cell(nC,1);
        Norma2 = cell(nC,1);
        nNeighbors = zeros(nC,1);
        p_values = zeros(nC,1);
        convergence_flags = false(nC,1);
        
        %% Main loop: barycentric regression network
        total_start_time = tic;
        
        for i = 1:nC
            fprintf('\n  Node %d / %d (%s) ---\n', i, nC, countries{i});
            
            y = pdfs2(:, i);
            Xall = pdfs2;
    
            
            % sort by similarity (smallest distance = most similar)
            [Wsort, pos] = sort(D(i, :));
            
            % remove self
            Wsort(pos == i) = [];
            pos(pos == i) = [];
            
            converged = false;
            
            % incrementally add neighbors
            for k = 1:min(maxK, numel(pos))
                X = Xall(:, pos(1:k));           % neighbor histograms
                S = size(X, 2);
                lambdav = repmat(1/S, S, 1);     % initial uniform weights
                
               % compute barycenter regression
                [P, lambdan, ~, Grad1, Grad2] = BarycenterRegression_Sinkhorn( ...
                    X, y, C, S, L, lambdav, maxiter, theta, alphav, epsilon, tol, tol2, dxmin);
                        
                P = P / sum(P);  % normalize
                
                % correct similarity test
                [pval, D_obs] = permtest_sinkhorn_divergence(y,P , C, epsilon, tol,maxiter,B);
                
                fprintf('    Using %2d neighbors → p = %.4f\n', k, pval);
                
                if pval >= alpha
                    fprintf('Converged at %d neighbors (p = %.4f)\n', k, pval);
                    AdjMatrix(pos(1:k), i) = lambdan;
                    Norma1{i} = Grad1;
                    Norma2{i} = Grad2;
                    nNeighbors(i) = k;
                    p_values(i) = pval;
                    convergence_flags(i) = true;
                    converged = true;
                    break;
                end
            end
            
            % if not converged, store last attempt anyway
            if ~converged
                fprintf(' No convergence for node %d (using %d neighbors)\n', i, k);
                AdjMatrix(pos(1:k), i) = lambdan;
                Norma1{i} = Grad1;
                Norma2{i} = Grad2;
                nNeighbors(i) = k;
                p_values(i) = pval;
                convergence_flags(i) = false;
            end
        end
        
        total_elapsed_time = toc(total_start_time);
        fprintf('\nTotal computation time: %.2f seconds\n', total_elapsed_time);
        
        %% Compute network properties
        % Create adjacency matrix from regression weights
        A = AdjMatrix > 0;  % Binary adjacency (any non-zero weight)
        
        % Number of nodes
        num_nodes = nC;
        
        % Number of links (directed)
        num_links = sum(A(:));
        
        % Network density (directed)
        max_possible_links = nC * (nC - 1);  % Directed, no self-loops
        density = num_links / max_possible_links;
        
        % Total flow (sum of all weights)
        total_flow = sum(abs(AdjMatrix(:)));
        
        % Average weight
        avg_weight = mean(AdjMatrix(AdjMatrix ~= 0));
        
        % In-degree and out-degree distributions
        out_degree = sum(A, 2);  % Number of outgoing connections
        in_degree = sum(A, 1)';  % Number of incoming connections
        
        avg_out_degree = mean(out_degree);
        avg_in_degree = mean(in_degree);
        
        % Convergence statistics
        convergence_rate = mean(convergence_flags) * 100;
        avg_neighbors = mean(nNeighbors);
        
        %% Compute additional network metrics
            
        % Average path length (using Dijkstra)
        avg_path_length = calculate_avg_path_length_directed(AdjMatrix);
        
        % Clustering coefficient (directed)
        clustering_coef = calculate_directed_clustering(AdjMatrix);
            
       
        
        %% Store results
        results(result_counter).alpha = alpha;
        results(result_counter).epsilon = epsilon;
        results(result_counter).AdjMatrix = AdjMatrix;
        results(result_counter).Norma1 = Norma1;
        results(result_counter).Norma2 = Norma2;
        results(result_counter).nNeighbors = nNeighbors;
        results(result_counter).p_values = p_values;
        results(result_counter).convergence_flags = convergence_flags;
        
        % Network properties
        results(result_counter).num_nodes = num_nodes;
        results(result_counter).num_links = num_links;
        results(result_counter).density = density;
        results(result_counter).total_flow = total_flow;
        results(result_counter).avg_weight = avg_weight;
        results(result_counter).avg_out_degree = avg_out_degree;
        results(result_counter).avg_in_degree = avg_in_degree;
        results(result_counter).convergence_rate = convergence_rate;
        results(result_counter).avg_neighbors = avg_neighbors;
        results(result_counter).avg_path_length = avg_path_length;
        results(result_counter).clustering_coefficient = clustering_coef;
        
        %% Display results for this parameter combination
        fprintf('\n  Network Properties:\n');
        fprintf('    Number of nodes: %d\n', num_nodes);
        fprintf('    Number of links: %.0f\n', num_links);
        fprintf('    Density: %.4f\n', density);
        fprintf('    Total flow: %.4f\n', total_flow);
        fprintf('    Average weight: %.4f\n', avg_weight);
        fprintf('    Average out-degree: %.2f\n', avg_out_degree);
        fprintf('    Average in-degree: %.2f\n', avg_in_degree);
        fprintf('    Convergence rate: %.1f%%\n', convergence_rate);
        fprintf('    Average neighbors: %.2f\n', avg_neighbors);
        fprintf('    Average path length: %.4f\n', avg_path_length);
        fprintf('    Clustering coefficient: %.4f\n', clustering_coef);
        
        %% Save individual results
        filename = sprintf('BarycentricNetwork_alpha%.3f_epsilon%.0e.mat', alpha, epsilon);
        save(filename, 'AdjMatrix', 'Norma1', 'Norma2', 'nNeighbors', ...
            'p_values', 'convergence_flags', 'results', 'countries');
        
        
        result_counter = result_counter + 1;
    end
end

%% Save all results
save('AllBarycentricResults.mat', 'results', 'countries', 'alpha_values', 'epsilon_values');

%% Create summary table
if ~isempty(results)
    summary_table = struct2table(results);
    disp('Summary of all regression network results:');
    disp(summary_table(:, {'alpha', 'epsilon', 'num_nodes', 'num_links', 'density', ...
        'total_flow', 'avg_weight', 'convergence_rate', 'avg_neighbors', ...
        'avg_path_length', 'clustering_coefficient'}));
end

%% ========================================================================
%% Helper function: Calculate average path length for directed graph
function avg_length = calculate_avg_path_length_directed(A)
% Calculate average shortest path length for directed weighted network
% A: weighted adjacency matrix (directed)

n = size(A, 1);
D = zeros(n);  % Distance matrix

% Convert weights to distances (higher weight = shorter distance)
% Use reciprocal for positive weights
D_temp = zeros(n);
for i = 1:n
    for j = 1:n
        if i ~= j
            if A(i,j) > 0
                D_temp(i,j) = 1 / A(i,j);
            else
                D_temp(i,j) = Inf;
            end
        else
            D_temp(i,j) = 0;
        end
    end
end

% Floyd-Warshall algorithm for directed graph
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
    for j = 1:n
        if i ~= j && D(i,j) < Inf
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
%% Helper function: Calculate directed clustering coefficient
function C = calculate_directed_clustering(A)
% Calculate clustering coefficient for directed graph
% Based on Fagiolo (2007)

n = size(A, 1);
C_i = zeros(n, 1);

% Create binary adjacency
A_bin = A > 0;

for i = 1:n
    % Find all neighbors (both successors and predecessors)
    successors = find(A_bin(i,:));
    predecessors = find(A_bin(:,i))';
    neighbors = unique([successors, predecessors]);
    k_i = length(neighbors);
    
    if k_i < 2
        C_i(i) = 0;
        continue;
    end
    
    % Count triangles
    triangles = 0;
    for idx1 = 1:k_i
        for idx2 = idx1+1:k_i
            n1 = neighbors(idx1);
            n2 = neighbors(idx2);
            
            % Check for any directed connection between n1 and n2
            if A_bin(n1,n2) || A_bin(n2,n1)
                triangles = triangles + 1;
            end
        end
    end
    
    % Possible triangles in directed graph
    possible_triangles = k_i * (k_i - 1);
    C_i(i) = triangles / possible_triangles;
end

C = mean(C_i);
end

%% ========================================================================
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

%% Barycentric Regression
function [P, lambdan, alphan, Norma1, Norma2] = BarycenterRegression_Sinkhorn( ...
        X, y, C, S, L, lambdav, maxiter, theta, alphav, gamma, tol, tol2, dxmin)
    %% BarycenterRegression_Sinkhorn
    %  Debiased Sinkhorn divergence loss with proper safeguards.

    %% Preallocate
    Norma1 = nan(maxiter, 1);
    Norma2 = nan(maxiter, 1);
    alphan = alphav;                                     % FIX 9: initialize properly

    N = size(C, 1);


    %% Initial run
    [~, wv] = MainOT_SinkhornDiv(X, y, N, C, S, L, lambdav, gamma);

    %% First update
    lambdan = projectSimplex(lambdav - alphav * wv);     % FIX 2: proper projection

    %% Second run
    [P, wn] = MainOT_SinkhornDiv(X, y, N, C, S, L, lambdan, gamma);

    %% Check convergence
    dx    = norm(lambdan - lambdav);
    gnorm = norm(wn - wv);
    
    if gnorm < tol || dx < dxmin
        alphan = alphav;
        return;                                          % FIX 9: explicit return
    end

    for t = 1:maxiter
        %% Adaptive step size
        grad_diff_norm = norm(wn - wv);
        alphan = min(sqrt(1 + theta) * alphav, ...
                         norm(lambdan - lambdav) / (2 * grad_diff_norm));
        

        lambdav = lambdan;
        wv      = wn;

        lambdan = projectSimplex(lambdav - alphan * wv); % FIX 2
        theta   = alphan / alphav;
        alphav  = alphan;                                % FIX 1: update alphav

        %% Compute barycenter and gradient
        [P, wn] = MainOT_SinkhornDiv(X, y, N, C, S, L, lambdan, gamma);

        %% Convergence diagnostics
        gnorm  = norm(wn - wv);
        gnorm2 = norm(wn);
        dx     = norm(lambdan - lambdav);

        Norma1(t) = gnorm;
        Norma2(t) = dx;

        if gnorm < tol || dx < dxmin
            break;
        end
        if gnorm2 < tol2
            break;
        end
    end
end


%% =========================================================================
%  MAIN OT: Barycenter + Sinkhorn Divergence gradient + Backward
%  =========================================================================
function [P, w] = MainOT_SinkhornDiv(X, y, N, C, S, L, lambda, gamma)

     %% Kernel (computed once, shared everywhere)
    K  = exp(-C / gamma);
    K(K < 1e-100) = 1e-100;
    KT = K';

    %% ----- Forward pass: Wasserstein barycenter P -----
    b   = ones(N, S, L);
    Phi = zeros(N, S, L);
    EPS = 1e-100;                                        % FIX 3: global floor

    for l = 2:L
        Phi_lambda = ones(N, 1);                         % FIX 4: start with ones
        for s = 1:S
            Kb = K * b(:, s, l-1);
            Kb = max(Kb, EPS);                           % FIX 3: prevent /0
            Phi(:, s, l)  = KT * (X(:, s) ./ Kb);
            Phi(:, s, l)  = max(Phi(:, s, l), EPS);     % FIX 3: prevent 0^lambda
            Phi_lambda    = Phi_lambda .* (Phi(:, s, l) .^ lambda(s));
        end
        P = Phi_lambda;
        b(:, :, l) = repmat(P, 1, S) ./ max(Phi(:, :, l), EPS);  % FIX 3
    end

    %% Normalize P to valid distribution
    P = max(P, EPS);
    P = P / sum(P);

    %% ----- Sinkhorn Divergence gradient w.r.t. P -----
    sinkhorn_tol  = 1e-10;
    sinkhorn_iter = 100;

    % (1) Cross term: OT_eps(P, y) — pass K to avoid recomputation  FIX 5
    [~, a_Py] = Sinkhorn_OT(C, gamma, P, y, sinkhorn_tol, sinkhorn_iter);

    % (2) Self term:  OT_eps(P, P)
    [~, a_PP] = Sinkhorn_OT(C, gamma, P, P, sinkhorn_tol, sinkhorn_iter);

    %% Gradient: dS/dP = eps*(log(a_Py) - log(a_PP)), times P for chain rule
    log_a_Py = log(max(a_Py, EPS));                      % FIX 6, 7: safe log
    log_a_PP = log(max(a_PP, EPS));

    g = gamma * (log_a_Py - log_a_PP) .* P;

    g(isnan(g)) = 0;
    g(isinf(g)) = 0;

    %% ----- Backward pass: propagate gradient to lambda -----
    w = zeros(S, 1);
    r = zeros(N, S);

    for l = L:-1:2
        for s = 1:S
            log_Phi = log(max(Phi(:, s, l), EPS));       % FIX 8: safe log
            w(s) = w(s) + log_Phi' * g;

            Phi_safe = max(Phi(:, s, l), EPS);           % FIX 8: prevent /0
            aus1 = (lambda(s) * g - r(:, s)) ./ Phi_safe;
            Kb   = max(K * b(:, s, l-1), EPS);           % FIX 8
            aus2 = X(:, s) ./ (Kb .^ 2);
            r(:, s) = -KT * (K * aus1 .* aus2) .* b(:, s, l-1);
        end
        g = sum(r, 2);
    end
end


%% =========================================================================
%  Sinkhorn_OT 
%  =========================================================================
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



%% =========================================================================
%  Simplex Projection (Duchi et al. 2008)
%  =========================================================================
function x = projectSimplex(v)
    %% Projects v onto {x >= 0, sum(x) = 1}
    n = length(v);
    u = sort(v, 'descend');
    cssv = cumsum(u);
    rho  = find(u .* (1:n)' > (cssv - 1), 1, 'last');
    if isempty(rho)                                      % safety
        x = ones(n, 1) / n;
        return;
    end
    tau = (cssv(rho) - 1) / rho;
    x   = max(v - tau, 0);
end



