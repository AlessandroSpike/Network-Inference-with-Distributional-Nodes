clc; clearvars; close all

%% ============================================================
%% Load data
load ReadyData.mat   % contains 'pdfs' and 'countries'

%% ============================================================
%% Parameters
nC      = numel(countries);
B       = 500;        % permutations
epsilon = 1e-5;
tol     = 1e-8;
niter   = 1000;
qFDR    = 0.01;       % FDR level

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
%% Initialize
P = ones(nC);
D = zeros(nC);
Errore = cell(nC);       % one error history per pair
%% ============================================================
%% MAIN LOOP (Sinkhorn divergence + permutation test)
for i = 1:nC
    for j = i+1:nC

        p1 = pdfs{i}(:);
        p2 = pdfs{j}(:);

        [pval, W_obs,Erro] = permtest_sinkhorn_divergence( ...
            p1, p2, C, epsilon, tol, niter, B);

        Errore{i,j}=Erro;
        P(i,j) = pval;
        P(j,i) = pval;

        D(i,j) = W_obs;
        D(j,i) = W_obs;

    end
end

%% ============================================================
%% Convert to similarity (RBF kernel)
sigma = median(D(D>0));
S = exp( -D.^2 ./ (2*sigma^2) );

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

%% ============================================================
%% Apply filtering
S_aus = S;
S_aus(P >= 0.01)= 0;
S_filtered = S;
S_filtered(P > p_cutoff) = 0;

%% ============================================================
%% Compare original vs FDR-filtered networks

% Ensure diagonals are zero
S_filtered(1:nC+1:end)      = 0;
S_aus(1:nC+1:end)  = 0;

%% ------------------------------------------------------------
%% Binary adjacency matrices
A_orig = S_aus > 0;
A_fdr  = S_filtered > 0;

% Upper triangle only (avoid double counting)
ut = triu(true(nC),1);

%% ------------------------------------------------------------
%% 1. Edge density
density_orig = sum(A_orig(ut)) / sum(ut(:));
density_fdr  = sum(A_fdr(ut))  / sum(ut(:));

%% ------------------------------------------------------------
%% 2. Total and average weight
total_weight_orig = sum(S_aus(ut));
total_weight_fdr  = sum(S_filtered(ut));

avg_weight_orig = mean(S_aus(A_orig & ut));
avg_weight_fdr  = mean(S_filtered(A_fdr & ut));

%% ------------------------------------------------------------
%% 3. Node strength comparison
strength_orig = sum(S_aus,2);
strength_fdr  = sum(S_filtered,2);

[r_strength,p_strength] = corr(strength_orig, strength_fdr, ...
                               'type','Spearman');

%% ------------------------------------------------------------
%% 4. Edge overlap (Jaccard index)
intersection = sum(A_orig & A_fdr & ut);
union_edges  = sum((A_orig | A_fdr) & ut);
jaccard_edge = intersection / union_edges;

%% ------------------------------------------------------------
%% 5. Weight correlation on common edges
common_edges = (A_orig & A_fdr & ut);

w_orig_common = S_aus(common_edges);
w_fdr_common  = S_filtered(common_edges);

if numel(w_orig_common) > 2
    [r_weight,p_weight] = corr(w_orig_common, w_fdr_common, ...
                               'type','Spearman');
else
    r_weight = NaN; p_weight = NaN;
end

%% ------------------------------------------------------------
%% 6. Frobenius distance between weighted adjacency matrices
fro_norm = norm(S_filtered - S_aus, 'fro');

%% ============================================================
%% Display results
fprintf('\n=== Network comparison ===\n');

fprintf('Edge density (original): %.3f\n', density_orig);
fprintf('Edge density (FDR):      %.3f\n\n', density_fdr);

fprintf('Total weight (original): %.3f\n', total_weight_orig);
fprintf('Total weight (FDR):      %.3f\n\n', total_weight_fdr);

fprintf('Avg weight (original):   %.3f\n', avg_weight_orig);
fprintf('Avg weight (FDR):        %.3f\n\n', avg_weight_fdr);

fprintf('Strength correlation (Spearman): r = %.3f, p = %.3g\n', ...
        r_strength, p_strength);

fprintf('Edge overlap (Jaccard):  %.3f\n', jaccard_edge);

fprintf('Weight corr. (common edges): r = %.3f, p = %.3g\n', ...
        r_weight, p_weight);

fprintf('Frobenius distance:      %.3f\n', fro_norm);

%% ============================================================
%% Save
save('SimilarityNetwork_FDR.mat', ...
    'S_filtered','S','S_aus','P','D','countries','Errore');

%% ============================================================
%% ============================================================
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
