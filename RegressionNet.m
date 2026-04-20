clc; clearvars; close all
load BarycentricRegressionNetwork.mat
load("ReadyData.mat",'countries')
%% Beautified Weighted Network Plot
G_sim = digraph(AdjMatrix, countries);
G_sim1 = digraph(AdjMatrix', countries);

figure('Color', [1 1 1], 'Position', [200 200 900 700]);

h = plot(G_sim, 'Layout','force');

% Node and edge aesthetics
h.MarkerSize = 7;
h.NodeColor = [0.2 0.6 0.8];
h.NodeLabelColor = [0 0 0];
h.LineWidth = 1.5;
h.EdgeAlpha = 0.7;
h.ArrowSize = 10;

% Scale edge widths by weights for better visibility
maxWidth = 3; % maximum line width
h.LineWidth = (G_sim.Edges.Weight / max(G_sim.Edges.Weight)) * maxWidth;

% Optionally, use a colormap for edges
colormap(jet)
edgeColors = G_sim.Edges.Weight;
h.EdgeCData = edgeColors;
colorbar

title('Weighted Regression Network', 'FontSize', 18, 'FontWeight', 'bold')

%% eigenvector
pg = centrality(G_sim,'pagerank');
pg2 = centrality(G_sim1,'pagerank');
indeg = centrality(G_sim,'indegree');
outdeg= centrality(G_sim,'outdegree');

figure
subplot(2,1,1)
bar([indeg,outdeg])
title('Degree')
grid on
axis tight
xticks(1:length(AdjMatrix))
xticklabels(countries)
legend('in-degree','out-degree')

subplot(2,1,2)
bar([pg,pg2])
title('PageRank')
grid on
axis tight
xticks(1:length(AdjMatrix))
xticklabels(countries)
legend('in-pagerank','out-pagerank')