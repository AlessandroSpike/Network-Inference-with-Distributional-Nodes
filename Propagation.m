clc;clearvars;close all
%% load
load("BarycentricRegressionNetwork.mat")
load ReadyData.mat   % contains 'pdfs' and 'countries'
nC       = numel(countries);
%% param
% entropy
epsilon=10^-5; 
% tolerance
tol=10^-5;
% num iterazioni
niter=50;
%% Define 2D cost matrix (bin positions)
edges = 0:1:10;
xgrid = edges(1:end-1) + 0.5;
ygrid = edges(1:end-1) + 0.5;
[X, Y] = meshgrid(xgrid, ygrid);
XY = [X(:), Y(:)];
C = pdist2(XY, XY, 'euclidean');
Ca = C / max(C(:));
%% centrality
centra =  centrality(digraph(AdjMatrix), 'pagerank','Importance',digraph(AdjMatrix).Edges.Weight);
% centra = sum(AdjMatrix);
%% Cascade Propagation with Level-by-Level Approach
threeD_Array = cat(3, pdfs{:});
pdfs2 = reshape(permute(threeD_Array, [3, 1, 2]), nC, [])';
pdfs2(pdfs2==0)=eps;
% Find initial shocked node
target_node = find(centra == min(centra));
target_node=target_node(1);
WP_afterShock = pdfs2;
WP_afterShock(:,target_node) = 1/size(pdfs2,1);
N = size(WP_afterShock, 2);

% Initialize tracking variables
shock = zeros(N, 1);          % Shock status (0 = no shock, 1 = shocked)
shock(target_node) = 1;       % Set initial shock
processed = zeros(N, 1);      % Track processed nodes
level_nodes = target_node;    % Nodes to process at current level
max_iterations = 200;         % Maximum iterations allowed

% Initialize dynamics tracking
num_shocked_per_step = zeros(max_iterations, 1);
num_shocked_per_step(1) = 1;  % Initial shocked node
shock_magnitude_per_step = zeros(max_iterations, 1);


[WassDist] = sinkhorn_divergence(WP_afterShock(:,target_node),pdfs2(:,target_node),Ca,epsilon, tol,niter);

shock_magnitude_per_step(1) = WassDist;


% Track unique shocked nodes at each step
unique_shocked_nodes_cumulative = zeros(max_iterations, 1);
unique_shocked_nodes_cumulative(1) = 1;  % Initial shocked node
currently_shocked = zeros(N, 1);
currently_shocked(target_node) = 1;

for step = 1:max_iterations
    if isempty(level_nodes)
        break;  % Stop if no more nodes to process
    end
    
    % Find all unprocessed neighbors of current level nodes
    next_level = [];
    shock_magnitude_this_step = 0;
    newly_shocked_this_step = zeros(N, 1);  % Track new shocks in this step
    
    for current_node = level_nodes'
        % Find valid neighbors
        valid_neighbors = [];
        neighbors_temp = find(AdjMatrix(current_node, :) > 0);
        
        % Only keep neighbors that are within bounds and unprocessed
        for n = 1:length(neighbors_temp)
            if neighbors_temp(n) <= N && processed(neighbors_temp(n)) == 0
                valid_neighbors = [valid_neighbors neighbors_temp(n)];
            end
        end
        
        % Process each valid neighbor
        for i = valid_neighbors
            % Skip if already shocked
            if currently_shocked(i) == 1
                continue;
            end
            
            % Find shocked neighbors that are connected to node i
            shocked_neighbors = [];
            temp_neighbors = find(AdjMatrix(i, :) > 0);
            
            % Check which neighbors are shocked
            for n = 1:length(temp_neighbors)
                if temp_neighbors(n) <= N && shock(temp_neighbors(n)) == 1
                    shocked_neighbors = [shocked_neighbors temp_neighbors(n)];
                end
            end
            
            if ~isempty(shocked_neighbors)
                % Store original distribution for magnitude calculation
                original_dist = WP_afterShock(:,i);
                
                % Prepare for Sinkhorn barycenter calculation
                DistrAttack_j = WP_afterShock(:, shocked_neighbors);
                C = cell(length(shocked_neighbors), 1);
                X = cell(length(shocked_neighbors), 1);
                
                link_w = AdjMatrix(i, shocked_neighbors);
                link_w=link_w/sum(link_w);
                
                
                % Setup cost matrices and distributions
                for c = 1:length(shocked_neighbors)
                    edges = 0:1:10;
                    xgrid = edges(1:end-1) + 0.5;
                    ygrid = edges(1:end-1) + 0.5;
                    [Xa, Ya] = meshgrid(xgrid, ygrid);
                    XY = [Xa(:), Ya(:)];
                    C2 = pdist2(XY, XY, 'euclidean');
                    C2 = C2 / max(C2(:));
                    C{c,1} = C2;
                    X{c,1} = DistrAttack_j(:,c);
                end
                
                % Calculate new distribution using Sinkhorn barycenter
                NewDistr = Sinkhorn_Barycenter(X, 100, length(C), 500, link_w, 1e-8);
                
                WP_afterShock(:,i) = NewDistr(:);
                
                % Calculate shock magnitude for this node
                [WassDist] = sinkhorn_divergence(NewDistr(:),original_dist,Ca,epsilon, tol,niter);
                shock_magnitude_this_step = shock_magnitude_this_step + WassDist;
                
                % Mark node as shocked and add to next level only if not already shocked
                if ~currently_shocked(i)
                    shock(i) = 1;
                    currently_shocked(i) = 1;
                    newly_shocked_this_step(i) = 1;
                    next_level = [next_level i];
                end
            end
        end
        processed(current_node) = 1;  % Mark current node as processed
    end
    
    % Update level_nodes for next iteration
    level_nodes = unique(next_level);
    
    % Record dynamics (only count newly shocked nodes)
    num_new_shocks = sum(newly_shocked_this_step);
    num_shocked_per_step(step+1) = num_new_shocks;
    unique_shocked_nodes_cumulative(step+1) = sum(currently_shocked);
    shock_magnitude_per_step(step+1) = shock_magnitude_this_step;
    
    % Break if all nodes are shocked
    if sum(currently_shocked) >= N
        break;
    end
end

% Trim zeros from tracking arrays
last_nonzero = find(num_shocked_per_step > 0, 1, 'last');
num_shocked_per_step = num_shocked_per_step(1:last_nonzero);
shock_magnitude_per_step = shock_magnitude_per_step(1:last_nonzero);
unique_shocked_nodes_cumulative = unique_shocked_nodes_cumulative(1:last_nonzero);

% Create visualization
figure


%Number of Nodes Shocked per Step
subplot(1,2,1)
yyaxis left
bar(1:length(num_shocked_per_step), num_shocked_per_step, 'FaceColor', [0.3 0.6 0.9])
ylabel('New Nodes Shocked', 'FontSize', 10)
yyaxis right
plot(1:length(unique_shocked_nodes_cumulative), unique_shocked_nodes_cumulative, 'r-', 'LineWidth', 2)
ylabel('Cumulative Shocked Nodes', 'FontSize', 10)
title('Propagation Dynamics: Shock Spread', 'FontSize', 12, 'FontWeight', 'bold')
xlabel('Step', 'FontSize', 10)
grid on
axis tight
legend('New Shocks', 'Cumulative', 'Location', 'northwest')
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
%Magnitude of Shock per Step
subplot(1,2,2)
yyaxis left
bar(1:length(shock_magnitude_per_step), shock_magnitude_per_step, 'FaceColor', [0.9 0.4 0.3])
ylabel('Magnitude of Change', 'FontSize', 10)
yyaxis right
plot(1:length(shock_magnitude_per_step), cumsum(shock_magnitude_per_step), 'b-', 'LineWidth', 2)
ylabel('Cumulative SChange', 'FontSize', 10)
title('Propagation Dynamics: Shock Magnitude per Step', 'FontSize', 12, 'FontWeight', 'bold')
xlabel('Step', 'FontSize', 10)
legend('Magnitude of Change', 'Cumulative', 'Location', 'northwest')
axis tight
grid on
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';






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


function P = Sinkhorn_Barycenter(X, N, S, L, lambda, gamma)
    % Computes the Wasserstein barycenter using the Sinkhorn-Knopp algorithm
    % 
    % Inputs:
    %   - X      : Cell array (Sx1) where each cell contains a (6x4) matrix
    %   - N      : Number of grid points (should be 24 for a 6x4 matrix)
    %   - S      : Number of input matrices (e.g., 59)
    %   - L      : Number of Sinkhorn iterations
    %   - lambda : Vector (Sx1) of weighting exponents for the barycenter computation
    %   - gamma  : Regularization parameter for Sinkhorn distance
    %
    % Output:
    %   - P : The Wasserstein barycenter (6x4 matrix)

    %% **Initialization**
    b = ones(N, S, L);      % Scaling matrix
    Phi = zeros(N, S, L);   % Transport matrix

    % Convert cell array to a 3D array (6x4xS)
    X = cat(3, X{:});  

    % Normalize each matrix to sum to 1 (convert to probability distributions)
    for s = 1:S
        X(:,:,s) = X(:,:,s) / sum(X(:,:,s), 'all');
    end

    % Reshape X into (24xS) for easier computation
    X = reshape(X, N, S);  
    X(X == 0) = 1e-50;      % Avoid numerical underflow

    %% **Construct Cost Matrix for 6x4 Field (24x24)**
    % Source points (e.g., grid in 2D)
    [X1, X2] = meshgrid(1:10, 1:10);     % 6x4 grid
    source_points = [X1(:), X2(:)];    % (24 x 2)
    
    % Target points (same or different grid)
    [Y1, Y2] = meshgrid(1:10, 1:10);     % Example target grid
    target_points = [Y1(:), Y2(:)];    % (24 x 2)
    C = pdist2(source_points, target_points, 'squaredeuclidean');
    C = C / max(C(:));  % Normalize the cost matrix

    %% **Precompute Kernels**
    K = exp(-C / gamma);
    K(K == 0) = 1e-50; % Avoid numerical issues
    KT = K';

    %% **Sinkhorn Iterations**
    for l = 2:L
        Phi_lambda = zeros(N, S);
        for s = 1:S
            % Compute transport updates
            Phi(:, s, l) = KT * (X(:, s) ./ (K * b(:, s, l-1)));
            Phi_lambda(:, s) = Phi(:, s, l).^lambda(s);
        end
        % Compute Wasserstein barycenter
        P = prod(Phi_lambda, 2);
        
        % Update scaling factors
        b(:,:,l) = repmat(P, 1, S) ./ Phi(:,:,l);
    end

    %% **Reshape Output**
    P = reshape(P, 10, 10); % Convert back to (6x4) matrix form
end



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
