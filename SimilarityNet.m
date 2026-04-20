clc; clearvars; close all
load SimilarityNetwork_FDR.mat
load("ReadyData.mat",'pdfs')

[Ci, Q] = modularity_und(S_filtered, 1);
clusters = unique(Ci);

%% Beautified Weighted Network Plot
G_sim = graph(S_filtered, countries);

figure('Color', [1 1 1], 'Position', [200 200 900 700]);

h = plot(G_sim, 'Layout', 'force');

% Node and edge aesthetics
h.MarkerSize = 7;
h.NodeColor = [0.2 0.6 0.8];
h.NodeLabelColor = [0 0 0];
h.LineWidth = 2;
h.EdgeAlpha = 0.7;

% Scale edge widths by weights for better visibility
maxWidth = 5; % maximum line width
h.LineWidth = (G_sim.Edges.Weight / max(G_sim.Edges.Weight)) * maxWidth;

title('Weighted Distance Network', 'FontSize', 18, 'FontWeight', 'bold')

% Optionally, use a colormap for edges
colormap(jet)
edgeColors = G_sim.Edges.Weight;
h.EdgeCData = edgeColors;
colorbar


%% Cluster distribution plots with balanced subplot layout
for i = 1:length(clusters)
    clusterIdx = find(Ci == clusters(i));
    nCountries = length(clusterIdx);
    
    % Determine balanced rows and columns
    nRows = ceil(sqrt(nCountries));
    nCols = ceil(nCountries / nRows);
    
    figure('Color', [1 1 1], 'Position', [150 150 1400 800]);
    
    for j = 1:nCountries
        subplot(nRows, nCols, j)
        imagesc(pdfs{clusterIdx(j)},'Interpolation','nearest');
        colormap('hot')
        colorbar
        axis tight
        title(countries{clusterIdx(j)}, 'FontSize', 12, 'FontWeight', 'bold')
        axis square
        xlabel('Immigration');
        ylabel('Economony');
    end
    
    sgtitle(['Distributions for Cluster ' num2str(clusters(i))], 'FontSize', 16, 'FontWeight', 'bold')
end


%% Parameters for Sinkhorn Barycenter
N = 100;                % 6x4 matrix flattened
L = 50;                % Number of Sinkhorn iterations
gamma = 1e-5;          % Regularization parameter

% Compute barycenters for each cluster
clusterBarycenters = cell(length(clusters),1);

for i = 1:length(clusters)
    clusterIdx = find(Ci == clusters(i));  % indices of countries in cluster
    S_cluster = length(clusterIdx);        % number of matrices in cluster
    
    if S_cluster == 0
        continue
    end
    
    % Extract matrices for the cluster
    X_cluster = pdfs(clusterIdx);
    
    % Equal weights for each matrix
    lambda = ones(S_cluster,1) / S_cluster;
    
    % Compute the barycenter using Sinkhorn_Barycenter
    P = Sinkhorn_Barycenter(X_cluster, N, S_cluster, L, lambda, gamma);
    
    clusterBarycenters{i} = P;
end



%% Visualize barycenters
figure('Color',[1 1 1], 'Position',[200 200 1000 400]);
xx = 1;
for i = length(clusters):-1:1
    if isempty(clusterBarycenters{i})
        continue
    end
    
    subplot(1,length(clusters),xx)
    h = bar3(clusterBarycenters{i});
     % --- This loop colors each bar according to its height (value) ---
    for k = 1:length(h)
        zdata = h(k).ZData;
        h(k).CData = zdata;
        h(k).FaceColor = 'interp';
    end
    
    % --- Apply original formatting ---
    colormap('hot');
    colorbar;
    axis tight;
    
    % Add labels for the 3D axes
    xlabel('Immigration Opinions');
    ylabel('Economony Opinions');
    zlabel('Probability');
    
    title(['Barycenter Cluster ' num2str(clusters(xx))], 'FontSize',14, 'FontWeight','bold');
    xx = xx+1;
end
sgtitle('Wasserstein Barycenters of Clusters','FontSize',16,'FontWeight','bold')


%% Function
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
    K(K == 0) = 1e-200; % Avoid numerical issues
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



function [Ci,Q]=modularity_und(A,gamma)
%MODULARITY_UND     Optimal community structure and modularity
%
%   Ci = modularity_und(W);
%   [Ci Q] = modularity_und(W,gamma);
%
%   The optimal community structure is a subdivision of the network into
%   nonoverlapping groups of nodes in a way that maximizes the number of
%   within-group edges, and minimizes the number of between-group edges.
%   The modularity is a statistic that quantifies the degree to which the
%   network may be subdivided into such clearly delineated groups.
%
%   Inputs:
%       W,
%           undirected weighted/binary connection matrix
%       gamma,
%           resolution parameter (optional)
%               gamma>1,        detects smaller modules
%               0<=gamma<1,     detects larger modules
%               gamma=1,        classic modularity (default)
%
%   Outputs:    
%       Ci,     optimal community structure
%       Q,      maximized modularity
%
%   Note:
%       This algorithm is essentially deterministic. The only potential
%       source of stochasticity occurs at the iterative finetuning step, in
%       the presence of non-unique optimal swaps. However, the present
%       implementation always makes the first available optimal swap and
%       is therefore deterministic.
%
%   References: 
%       Newman (2006) -- Phys Rev E 74:036104, PNAS 23:8577-8582.
%       Reichardt and Bornholdt (2006) Phys Rev E 74:016110.
%
%   2008-2016
%   Mika Rubinov, UNSW
%   Jonathan Power, WUSTL
%   Dani Bassett, UCSB
%   Xindi Wang, Beijing Normal University
%   Roan LaPlante, Martinos Center for Biomedical Imaging

%   Modification History:
%   Jul 2008: Original (Mika Rubinov)
%   Oct 2008: Positive eigenvalues made insufficient for division (Jonathan Power)
%   Dec 2008: Fine-tuning made consistent with Newman's description (Jonathan Power)
%   Dec 2008: Fine-tuning vectorized (Mika Rubinov)
%   Sep 2010: Node identities permuted (Dani Bassett)
%   Dec 2013: Gamma resolution parameter included (Mika Rubinov)
%   Dec 2013: Detection of maximum real part of eigenvalues enforced (Mika Rubinov)
%               Thanks to Mason Porter and Jack Setford, University of Oxford
%   Dec 2015: Single moves during fine-tuning enforced (Xindi Wang)
%   Jan 2017: Removed node permutation and updated documentation (Roan LaPlante)

if ~exist('gamma','var')
    gamma = 1;
end

N=length(A);                            %number of vertices
% n_perm = randperm(N);                   %DB: randomly permute order of nodes
% A = A(n_perm,n_perm);                   %DB: use permuted matrix for subsequent analysis
K=sum(A);                               %degree
m=sum(K);                               %number of edges (each undirected edge is counted twice)
B=A-gamma*(K.'*K)/m;                    %modularity matrix
Ci=ones(N,1);                           %community indices
cn=1;                                   %number of communities
U=[1 0];                                %array of unexamined communites

ind=1:N;
Bg=B;
Ng=N;

while U(1)                              %examine community U(1)
    [V,D]=eig(Bg);
    [~,i1]=max(real(diag(D)));          %maximal positive (real part of) eigenvalue of Bg
    v1=V(:,i1);                         %corresponding eigenvector

    S=ones(Ng,1);
    S(v1<0)=-1;
    q=S.'*Bg*S;                         %contribution to modularity

    if q>1e-10                       	%contribution positive: U(1) is divisible
        qmax=q;                         %maximal contribution to modularity
        Bg(logical(eye(Ng)))=0;      	%Bg is modified, to enable fine-tuning
        indg=ones(Ng,1);                %array of unmoved indices
        Sit=S;
        while any(indg)                 %iterative fine-tuning
            Qit=qmax-4*Sit.*(Bg*Sit); 	%this line is equivalent to:
            [qmax,imax]=max(Qit.*indg); %for i=1:Ng
            Sit(imax)=-Sit(imax);       %	Sit(i)=-Sit(i);
            indg(imax)=nan;             %	Qit(i)=Sit.'*Bg*Sit;
            if qmax>q                   %	Sit(i)=-Sit(i);
                q=qmax;                 %end
                S=Sit;
            end
        end

        if abs(sum(S))==Ng              %unsuccessful splitting of U(1)
            U(1)=[];
        else
            cn=cn+1;
            Ci(ind(S==1))=U(1);         %split old U(1) into new U(1) and into cn
            Ci(ind(S==-1))=cn;
            U=[cn U];                   %#ok<AGROW>
        end
    else                                %contribution nonpositive: U(1) is indivisible
        U(1)=[];
    end

    ind=find(Ci==U(1));                 %indices of unexamined community U(1)
    bg=B(ind,ind);
    Bg=bg-diag(sum(bg));                %modularity matrix for U(1)
    Ng=length(ind);                     %number of vertices in U(1)
end

s=Ci(:,ones(1,N));                      %compute modularity
Q=~(s-s.').*B/m;
Q=sum(Q(:));
% Ci_corrected = zeros(N,1);              % DB: initialize Ci_corrected
% Ci_corrected(n_perm) = Ci;              % DB: return order of nodes to the order used at the input stage.
% Ci = Ci_corrected;                      % DB: output corrected community assignments
end


