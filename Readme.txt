Project: Cross-Country Socio-Economic Analysis

This project processes and analyzes public opinion data from the European Social Survey (ESS11). It focuses on understanding the relationship between public satisfaction with the economy and perceptions of immigration across different countries.

===============================================================================
Script 1: DataPreprocess.m
===============================================================================

Purpose
-------
The primary goal of this script is to clean raw survey data and transform individual-level responses into 2D probability distributions (joint PDFs) for each country. This allows for a comparative analysis of how two specific variables correlate within different national contexts.

Data Source
-----------
- Input File: ESS11.csv (European Social Survey Round 11).
- Variables Extracted:
  - cntry: The country code.
  - stfeco: Satisfaction with the state of the economy (Scale 0-10).
  - imbgeco: Perception of whether immigration is good or bad for the economy (Scale 0-10).

Workflow
--------
1. Data Cleaning:
   - Reads the CSV file into a table.
   - Filters out missing values (NaN), invalid country entries, and responses outside the 0-10 range.
2. Grid Generation:
   - Defines a 10x10 grid (using edges from 0 to 10) to represent the possible combinations of responses.
3. Distribution Calculation (PDF):
   - Iterates through each unique country in the dataset.
   - Uses histcounts2 to calculate a 2D histogram of the two variables.
   - Normalization: Converts raw counts into a probability distribution where the sum of all bins equals 1.
4. Storage:
   - Saves the processed distributions and the list of country names for use in subsequent analysis scripts.

Requirements
------------
- ESS11.csv located in the same directory.

Output
------
- ReadyData.mat: A MATLAB workspace file containing:
  - pdfs: A cell array where each cell contains a 10x10 matrix representing the joint probability distribution for a specific country.
  - countries: A list of unique country codes corresponding to the pdfs.

===============================================================================
Script 2: Preliminary.m
===============================================================================

Purpose
-------
This script generates professional, sorted, diverging stacked bar charts. It breaks down the 2D joint distributions into two separate 1D marginal analyses:
1. Immigration Opinions: How respondents feel about immigration's impact on the economy.
2. Economic Satisfaction: How satisfied respondents are with the current state of their national economy.

Key Features
------------
- Marginalization: The script uses the sum function to "integrate out" one variable, allowing the analysis of immigration and the economy independently.
- Diverging Stacked Bar Charts:
  - Negative Sentiment (Scale 0-4): Plotted to the left of the zero axis (Blue/Navy tones).
  - Positive Sentiment (Scale 5-10): Plotted to the right of the zero axis (Orange/Red tones).
- Dynamic Sorting: Countries are automatically sorted by the highest percentage of "Positive" responses, making it easy to identify the most optimistic vs. pessimistic nations.
- Custom Aesthetics: Uses a specific color palette sampled for high readability and removes unnecessary chart junk (like black edges and boxes).

Visualization Logic
-------------------
- The 11-point scale (0-10) is split into two groups:
  - Indices 1-4: Grouped as Negative/Conservative sentiment.
  - Indices 5-10: Grouped as Positive/Liberal sentiment.
- The X-axis is normalized to show percentages (%), centered at 0 to highlight the balance of opinion.

How to Run
----------
1. Ensure ESS11.csv is in the project folder.
2. Run DataPreprocess.m to generate the ReadyData.mat file.
3. Run Preliminary.m to generate the two comparison figures:
   - Figure 1: Immigration Opinions by Country.
   - Figure 2: Economic Opinions by Country.

Requirements
------------
- Statistics and Machine Learning Toolbox (for histcounts2)

===============================================================================
Script 3: MainSimilarity.m
===============================================================================

Purpose
-------
This script identifies which countries are statistically "similar" in their 2D opinion profiles. Instead of just comparing means, it compares the entire shape of the distributions using Optimal Transport.

Mathematical Background
-----------------------
- Cost Matrix (C): Defines the "distance" between bins on the 0-10 scale. A person answering "4" is mathematically closer to "5" than to "10."
- Sinkhorn Divergence: An efficient version of the Earth Mover's Distance. It measures the "work" required to transform one country's opinion distribution into another's.
- Permutation Test: For every pair of countries (i, j), the script performs B=500 random shuffles to determine if the observed similarity is statistically significant or happened by chance.
- RBF Kernel: Converts the calculated distances into a Similarity Score (S) ranging from 0 to 1.

Network Filtering & FDR
-----------------------
Comparing many countries leads to the "Multiple Comparisons Problem." The script implements:
1. Benjamini-Hochberg (FDR) Procedure: A rigorous statistical thresholding method to control the number of false-positive similarities.
2. Network Comparison: The script compares the "Original" network (all correlations) against the "Filtered" network (only statistically significant connections).

Key Metrics Calculated
----------------------
- Edge Density: How connected the network of countries is.
- Node Strength: Which countries are "hubs" of typical/average opinion.
- Jaccard Index: The overlap between the original and FDR-filtered network structures.
- Frobenius Distance: The mathematical difference between the weighted adjacency matrices of the two networks.

Requirements
------------
- Parallel Computing Toolbox: Used for parfor to speed up the 500-iteration permutation tests.
- ReadyData.mat: Must be generated by DataPreprocess.m.

Output
------
- SimilarityNetwork_FDR.mat: Contains the full similarity matrix (S), the filtered matrix (S_filtered), and p-values (P) for all country pairings.
- Console Log: Detailed comparison of network density, weight correlations, and FDR cutoffs.

===============================================================================
Script 4: PlotErrorSimil.m
===============================================================================

Purpose
-------
When calculating the similarity between two countries, the script uses an iterative algorithm (Sinkhorn-Knopp). This script acts as a Quality Assurance step to ensure that the algorithm converged correctly across all country-pair comparisons.

Key Features
------------
- Data Aggregation: Extracts the error history from the Errore cell array (generated in MainSimilarity.m).
- Statistical Summary: Calculates the Mean Error and Standard Deviation across all hundreds of country-pair optimizations.
- Convergence Visualization:
  - Plots the average error decrease over iterations.
  - Includes a shaded region representing the standard deviation (variability) of the convergence.

Why This Is Important
---------------------
If the error curve does not flatten out or remain low, it suggests that the "Similarity" scores calculated in the project might be numerically unstable. A smooth, downward-sloping curve confirms that the results are mathematically sound.

Requirements
------------
- SimilarityNetwork_FDR.mat: Must be generated by MainSimilarity.m.

===============================================================================
Script 5: SimilarityNet.m
===============================================================================

Purpose
-------
This script analyzes the filtered country similarity network, detects clusters of similar countries, and visualizes both the network structure and the representative distributions of each cluster.

Workflow
--------
1. Loads SimilarityNetwork_FDR.mat and ReadyData.mat.
2. Applies modularity-based community detection to the filtered similarity matrix S_filtered.
3. Builds and plots a weighted country network using a force-directed layout.
4. Groups countries by cluster and displays their 10x10 joint PDFs with imagesc.
5. Computes a Wasserstein barycenter for each cluster using a custom Sinkhorn algorithm.
6. Visualizes each cluster barycenter as a 3D bar plot.

Requirements
------------
- SimilarityNetwork_FDR.mat
- ReadyData.mat

Output
------
- Weighted similarity network figure
- One figure per detected cluster showing country distributions
- Final figure with Wasserstein barycenters for each cluster

===============================================================================
Script 6: SimilarityThresoldComparison.m
===============================================================================

Purpose
-------
This script compares different network pruning strategies applied to the country similarity matrix, showing how each thresholding method changes the final graph structure.

Workflow
--------
1. Loads SimilarityNetwork_FDR.mat and sets the significance level alpha = 0.01.
2. Rebuilds three statistically filtered networks: Original (P < alpha), FDR (Benjamini-Hochberg), and FWER (Bonferroni).
3. Creates additional quantile-based pruned networks by removing edges above selected distance thresholds.
4. Converts all weighted networks into binary adjacency matrices for structural comparison.
5. Computes graph metrics for each method: edge density, total weight, average weight, node strength, Jaccard overlap, weight correlation, and Frobenius distance.
6. Prints a full comparison report and summary table in the MATLAB console.
7. Produces bar plots comparing density, overlap with the original network, and average retained edge weight.

Requirements
------------
- SimilarityNetwork_FDR.mat

Output
------
- Console summary of threshold values and network comparison metrics
- Figure with three bar charts comparing pruning methods

===============================================================================
Script 7: MainRegression.m
===============================================================================

Purpose
-------
This script builds a barycentric regression network between countries by expressing each country's distribution as a Wasserstein barycenter of its most similar neighbors.

Workflow
--------
1. Loads ReadyData.mat and the precomputed distance matrix D from SimilarityNetwork_FDR.mat.
2. Builds a 2D transport cost matrix for the 10x10 opinion grid and reshapes all country PDFs into column vectors.
3. For each country, sorts all other countries by increasing distance and incrementally adds neighbors.
4. Runs Sinkhorn-based barycenter regression to estimate the optimal neighbor weights.
5. Uses a permutation test on the Sinkhorn divergence to decide when the reconstructed barycenter is statistically similar to the target country.
6. Stores the regression weights, convergence diagnostics, and number of neighbors used for each node.
7. Saves the final results in BarycentricRegressionNetwork.mat.

Requirements
------------
- ReadyData.mat
- SimilarityNetwork_FDR.mat
- Parallel Computing Toolbox recommended (parfor used in permutation tests)

Output
------
- BarycentricRegressionNetwork.mat containing AdjMatrix, Norma1, Norma2, and nNeighbors

===============================================================================
Script 8: PlotErrorRegress.m
===============================================================================

Purpose
-------
This script evaluates the convergence behavior of the barycentric regression algorithm by averaging the error histories stored during optimization.

Workflow
--------
1. Loads Norma1 and Norma2 from BarycentricRegressionNetwork.mat.
2. Extracts all non-empty error curves from Norma1, which track gradient norm changes across iterations.
3. Computes the mean and standard deviation of these curves over all country regressions.
4. Plots the average gradient error with a shaded standard deviation band.
5. Repeats the same procedure for Norma2, which tracks changes in the regression coefficients.
6. Produces a second subplot showing the average coefficient-update error with variability.

Requirements
------------
- BarycentricRegressionNetwork.mat

Output
------
- One figure with two panels:
  - Left: Average gradient norm decay
  - Right: Average regression coefficient norm decay

===============================================================================
Script 9: RegressionNet.m
===============================================================================

Purpose
-------
This script visualizes the directed barycentric regression network and summarizes the centrality of each country within that network.

Workflow
--------
1. Loads BarycentricRegressionNetwork.mat and the country labels from ReadyData.mat.
2. Builds a directed graph from AdjMatrix, where edge weights represent regression coefficients.
3. Plots the weighted directed network with edge widths and colors scaled by coefficient strength.
4. Computes node centrality measures, including in-degree, out-degree, and PageRank.
5. Repeats PageRank on the transposed network to compare incoming vs outgoing importance.
6. Displays bar charts comparing degree and PageRank values across countries.

Requirements
------------
- BarycentricRegressionNetwork.mat
- ReadyData.mat

Output
------
- Weighted directed regression network figure
- Centrality figure with degree and PageRank comparisons

===============================================================================
Script 10: Propagation.m
===============================================================================

Purpose
-------
This script simulates how a shock spreads through the directed barycentric regression network, starting from the most central country and propagating distributional changes level by level.

Workflow
--------
1. Loads BarycentricRegressionNetwork.mat and the country PDFs from ReadyData.mat.
2. Computes a 2D transport cost matrix and PageRank centrality on the regression network.
3. Selects the most central node as the initial shocked country and replaces its distribution with a uniform one.
4. Propagates the shock through outgoing network links, updating each affected country via a weighted Sinkhorn barycenter of shocked neighbors.
5. Measures the size of each update using Sinkhorn divergence between old and new distributions.
6. Tracks the number of newly shocked nodes, cumulative shocked nodes, and total shock magnitude at each propagation step.
7. Visualizes the spread of the cascade and the cumulative magnitude of change over time.

Requirements
------------
- BarycentricRegressionNetwork.mat
- ReadyData.mat

Output
------
- Figure with two panels:
  - Left: New and cumulative shocked nodes per step
  - Right: Shock magnitude per step and cumulative change

===============================================================================
Script 11: MainSimilarityRobustness.m
===============================================================================

Purpose
-------
This script performs a robustness analysis of the similarity network by systematically varying two critical algorithmic parameters: the statistical significance threshold (alpha) and the Sinkhorn regularization strength (epsilon). It evaluates how network topology and statistical properties change across parameter combinations.

Workflow
--------
1. Parameter Grid Search: Tests all combinations of 3 significance levels (alpha = 0.01, 0.05, 0.1) and 3 regularization values (epsilon = 1e-4, 1e-5, 1e-6).
2. Statistical Testing: For each configuration, recalculates country-pair similarities using permutation tests with Sinkhorn divergence.
3. FDR Filtering: Applies Benjamini-Hochberg correction to control false discovery rate at each alpha level.
4. Network Metrics: Computes density, total flow, average path length, and weighted clustering coefficient for each resulting network.
5. Comparative Output: Generates a summary table showing how network structure varies with parameter choices.

Requirements
------------
- ReadyData.mat
- Parallel Computing Toolbox (for parfor in permutation tests)

Output
------
- Individual Files: SimilarityNetwork_alpha[X]_epsilon[Y].mat for each parameter combination (9 files total)
- AllResults.mat: Consolidated results structure with all network metrics
- Console Table: Summary comparing density, clustering, and path length across all 9 configurations

===============================================================================
Script 12: MainRegressionRobustness.m
===============================================================================

Purpose
-------
This script extends the barycentric regression framework by testing its stability across different statistical thresholds and regularization parameters. It evaluates how the directed regression network structure responds to variations in significance level (alpha) and Sinkhorn entropy (epsilon).

Workflow
--------
1. Parameter Sweep: Tests all combinations of 3 significance thresholds (alpha = 0.01, 0.05, 0.1) and 3 regularization values (epsilon = 1e-4, 1e-5, 1e-6).
2. Adaptive Neighbor Selection: For each country and parameter set, incrementally adds neighbors (sorted by similarity) until the barycentric reconstruction passes a permutation test at the given alpha level.
3. Convergence Tracking: Records gradient norms, coefficient updates, and the number of neighbors required for statistical equivalence.
4. Network Analysis: Computes directed graph metrics including in/out-degree distributions, average path length, and directed clustering coefficient.
5. Comparative Output: Generates a summary table showing convergence rates, network density, and topological properties across all 9 configurations.

Requirements
------------
- ReadyData.mat
- SimilarityNetwork_alpha0.010_epsilon[X].mat files (precomputed distance matrices from MainSimilarityRobustness.m)
- Parallel Computing Toolbox (for permutation tests)

Output
------
- Individual Files: BarycentricNetwork_alpha[X]_epsilon[Y].mat for each parameter combination (9 files total)
- AllBarycentricResults.mat: Consolidated structure with adjacency matrices, convergence diagnostics, and network metrics
- Console Table: Comparison of convergence success rates, average neighbors used, and directed clustering

===============================================================================
Script 13: SimilarityThresholdComparison.m
===============================================================================

Purpose
-------
This script performs a systematic comparison of network pruning strategies applied to the country similarity matrix. It evaluates how different statistical correction methods and distance-based thresholds affect the resulting network structure, providing insight into the robustness of identified country relationships.

Workflow
--------
1. Statistical Thresholding: Constructs three networks using different multiple comparison corrections:
   - Original: Raw p-value threshold (alpha = 0.01, no correction)
   - FDR: Benjamini-Hochberg procedure (controls false discovery rate)
   - FWER: Bonferroni correction (controls family-wise error rate)
2. Quantile-Based Pruning: Creates five additional networks by removing edges above distance quantiles (1%, 5%, 10%, 25%, 50%), independent of statistical significance.
3. Structural Comparison: For all 8 pruning methods, computes:
   - Edge density and total network weight
   - Jaccard overlap with the original network
   - Spearman correlation of node strength and edge weights
   - Frobenius distance between adjacency matrices
4. Visualization: Generates bar charts comparing density, overlap, and average edge weight across all pruning strategies.

Requirements
------------
- SimilarityNetwork_FDR.mat

Output
------
- Console Report: Detailed table showing edge counts, density, Jaccard indices, and correlations for all 8 methods
- Figure: Three-panel comparison of network properties (density, overlap, average weight)

===============================================================================
Script 14: KL_comparison.m
===============================================================================

Purpose
-------
This script provides a rigorous methodological justification for using Wasserstein distance instead of Kullback-Leibler (KL) divergence in the similarity analysis. It addresses three fundamental questions through controlled experiments: metric validity, geometric sensitivity, and network structural differences.

Workflow
--------
Experiment 1: Metric Properties
- Tests the triangle inequality on 50 randomly selected country triplets.
- Demonstrates that KL divergence violates the mathematical properties of a proper distance metric, while Wasserstein distance does not.

Experiment 2: Geometric Sensitivity
- Creates two synthetic transformations of a country's opinion distribution:
  - Spatial Shift: Moves opinion mass systematically across the grid (e.g., rightward shift in economic satisfaction).
  - Random Permutation: Scrambles the same mass without preserving spatial structure.
- Compares how each metric responds to meaningful vs. meaningless changes.

Experiment 3: Network Structure Comparison
- Computes full pairwise distance matrices using both Wasserstein and symmetrized KL divergence.
- Calculates Spearman and Pearson correlations between the two distance measures.
- Builds similarity networks using the top 20% of edges from each method.
- Identifies country pairs where the two metrics disagree most in ranking.

Requirements
------------
- ReadyData.mat
- MATLAB Statistics and Machine Learning Toolbox

Output
------
- Wasserstein_vs_KL_FocusedComparison.mat: Distance matrices, adjacency matrices, and comparison metrics
- Console Report: Triangle inequality violations, sensitivity ratios, correlation statistics, and top rank disagreements
- Figure: Three-panel visualization showing distance correlation, sensitivity comparison, and network edge overlap (Venn diagram)

===============================================================================
Script 15: RBF_comparison.m
===============================================================================

Purpose
-------
This script provides a comprehensive comparison of distance metrics and similarity transformations used in the network analysis. It systematically evaluates three alternative approaches: (1) Wasserstein distance with RBF kernel, (2) Euclidean distance with RBF kernel, and (3) raw Wasserstein distance with inverse transformation, addressing methodological choices in converting distances to similarity weights.

Workflow
--------
Distance Computation
- Calculates pairwise Wasserstein distances (Sinkhorn divergence) between all country distributions.
- Calculates pairwise Euclidean distances between probability vectors.

Similarity Transformations
- Converts distances to similarities using three methods:
  - Wasserstein + RBF: exp(-d^2 / (2*sigma^2)) with median-based bandwidth selection
  - Euclidean + RBF: Same kernel applied to Euclidean distances
  - Wasserstein Raw: 1/(1+d) transformation without kernelization

Comparative Analysis
- Generates nine detailed tables:
  1. Distance Statistics: Mean, median, range, and correlation between Wasserstein and Euclidean
  2. Similarity Statistics: Distribution properties of all three similarity measures
  3. Pairwise Correlations: Spearman correlations between similarity methods
  4. Node Strength: Total connectivity per country under each method
  5. Top Countries: Ranking of most central countries (highest node strength)
  6. Discordant Pairs: Country pairs where Wasserstein and Euclidean disagree most
  7. Concordant Pairs: Pairs where both metrics agree strongly
  8. Similarity Range: Effect of RBF kernel on dynamic range
  9. Percentile Analysis: Distribution of distances and similarities at key quantiles

Requirements
------------
- ReadyData.mat
- MATLAB Statistics and Machine Learning Toolbox

Output
------
- DistanceMetric_Comparison.mat: All three similarity matrices, distance matrices, and statistical comparison tables
- Console Report: Nine formatted tables with detailed statistical comparisons
- Figure 1: Distance vs. similarity correlations between Wasserstein and Euclidean methods
- Figure 2: Comparison of raw inverse transformation vs. RBF kernelization