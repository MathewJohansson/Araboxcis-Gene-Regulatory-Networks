



# Project: Gene Regulatory Network Analysis
# Focus: G-box regulatory networks in flower development of Arabidopsis thaliana


# Required packages
# install.packages('Matrix')
# install.packages('umap')
# install.packages('pheatmap')
# install.packages("BiocManager")
library(Matrix)
library(umap)
library(pheatmap)

# Load helper functions
source('dev/utilities/dataprocessingHelperFunctions.R')

# Create output directories if they don't exist
dir.create("results", showWarnings = FALSE)
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/networks", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)




# 1. LOAD DATA -----------------------------------------------------------------

# Load the flower single-cell RNA-seq dataset
load("data/flowerdata.RData", verbose = TRUE)
a <- load("data/flowerdata.RData", verbose = TRUE)

# Load the original AraBOXcis network (trained on bulk RNA-seq in seedlings)
araboxcis <- read.csv(file = 'data/gboxNetwork22C.csv', header = TRUE)




# 2. EXPLORE DATASET CONTENTS --------------------------------------------------

# Check what objects were loaded
print(a)

# Dimensions of gbox (gene expression matrix)
dim(gbox)

# Number of clusters
length(clust)

# Gene names (first 10)
rownames(gbox)[1:10]

# Cell names (first 10)
colnames(gbox)[1:10]

# Preview gene expression values
as.matrix(gbox[1:50, 1:3])




# 3. EXAMINE CELL CLUSTERS -----------------------------------------------------

# Cluster designations of the first 8 cells
clust[1:8]

# Cell counts in each cluster
table(clust)

# Visualize cell distribution across clusters
plot(table(sort(clust)), 
     xlab = 'Cluster Name', 
     ylab = 'Number of Cells', 
     main = 'Cell Distribution Across Clusters - Flower Data')




# 4. EXAMINE ARABOXCIS NETWORK -------------------------------------------------

# Dimensions of the AraBOXcis network
dim(araboxcis)

# First 4 edges of the network
araboxcis[1:4, ]

# Distribution of edge scores
hist(araboxcis[, 3],
     main = 'Distribution of AraBOXcis Edge Scores',
     xlab = 'Edge Score')

# Extract transcription factors from AraBOXcis network
tfs <- unique(araboxcis[, 1])

# Filter for TFs present in the single-cell dataset
tfSubs <- tfs[which(tfs %in% rownames(gbox))]

# Number of TFs in both datasets
length(tfSubs)




# 5. FILTER GENES AND CELLS WITH LOW EXPRESSION --------------------------------

# Original dimensions
dim(gbox)

# Filter out cells with <1% of genes expressed
thresh <- 0.01
numberGenesPerCell <- apply(gbox, 
                            2, 
                            function(i) {
                              length(which(i > 0))
                            })

includeCells <- which(numberGenesPerCell > (thresh * dim(gbox)[1])) 
gbox_filtered <- gbox[, includeCells] 

cat("After cell filtering:", dim(gbox_filtered)[2], "cells remaining\n")

# Remove genes expressed in <1% of cells
numberCellPerGene <- apply(gbox, 
                           1, 
                           function(i) {
                             length(which(i > 0))
                           })

includeGenes <- which(numberCellPerGene > (thresh * dim(gbox)[2]))
gbox_filtered <- gbox_filtered[includeGenes, ]

cat("After gene filtering:", dim(gbox_filtered)[1], "genes remaining\n")




# 6. DIMENSIONALITY REDUCTION - UMAP -------------------------------------------

# UMAP on filtered data
gbox_umap <- umap(gbox_filtered)

# Visualize with cluster colors
colours <- rainbow(length(unique(clust)))

plot(gbox_umap$layout[, 1], 
     gbox_umap$layout[, 2],
     col = colours[clust[includeCells]], 
     pch = 20,
     main = 'UMAP - Flower Cell Clusters', 
     xlab = 'UMAP Component 1', 
     ylab = 'UMAP Component 2')




# 7. PCA FOLLOWED BY UMAP ------------------------------------------------------

# Perform PCA (retain 5 components)
pca <- prcomp(gbox_filtered, 
              scale. = TRUE, 
              rank. = 5)

# UMAP on PCA results
gbox_pca_umap <- umap(pca$x)

plot(gbox_pca_umap$layout[, 1], 
     gbox_pca_umap$layout[, 2], 
     col = colours[clust[includeCells]], 
     pch = 20, 
     main = 'PCA + UMAP - Flower Cell Clusters', 
     xlab = 'UMAP Component 1', 
     ylab = 'UMAP Component 2')




# 8. GENE EXPRESSION WITHIN CLUSTERS -------------------------------------------

# Convert cluster assignments to numeric
clustAsNumbers <- as.numeric(paste(clust))

# Calculate average expression of each gene in each cluster
geneExpByCluster <- apply(gbox,
                          1,
                          function(i) {
                            sapply(0:(length(unique(clust)) - 1), 
                                   function(j) {
                                     ids <- which(clustAsNumbers == j)
                                     mean(i[ids])
                                   })
                          })

colnames(geneExpByCluster) <- rownames(gbox)

# Dimensions: rows = clusters, columns = genes
dim(geneExpByCluster) 




# 9. VISUALIZE GENE EXPRESSION AS HEATMAP --------------------------------------

# Heatmap of average gene expression per cluster (scaled by column)
pheatmap(geneExpByCluster, 
         scale = 'column',
         main = 'Average Gene Expression by Cluster')

# Load cluster labels (cell type annotations)
clustLabs <- read.table('data/clusterLabels.txt', 
                        header = TRUE, 
                        sep = '\t')

# Check available organs in the dataset
unique(clustLabs[, 'Organ'])

# Extract flower-specific cluster names
organ <- 'Flower'
simpleNames <- clustLabs[which(clustLabs[, 'Organ'] == organ), "Cell.type.suggested"]
print(simpleNames)

# Set meaningful row names for the heatmap
rownames(geneExpByCluster) <- simpleNames

# Redraw heatmap with cell type labels
pheatmap(geneExpByCluster, 
         scale = 'column',
         main = 'Average Gene Expression by Cell Type - Flower')

# Save the table for future use
write.table(geneExpByCluster, 
            'results/tables/Flower_avgExpressionByCluster.txt',
            sep = '\t',
            quote = FALSE)




# 10. FOCUS ON TRANSCRIPTION FACTORS -------------------------------------------

# Heatmap focusing only on G-box binding TFs
pheatmap(geneExpByCluster[, tfSubs], 
         scale = 'column',
         main = 'TF Expression Across Cell Types')




# 11. CORRELATION BETWEEN TF AND TARGET GENE EXPRESSION ------------------------

# Calculate Spearman correlation between every TF and every gene
# This examines co-expression patterns at the cell type level
corMat <- sapply(tfSubs, function(tf) {
  apply(geneExpByCluster, 
        2,
        function(gene) {
          cor(geneExpByCluster[, tf], gene, method = 'spearman')
        })
})

# Dimensions: rows = genes, columns = TFs
dim(corMat)

# Visualize TF-gene correlation patterns
pheatmap(corMat,
         main = 'TF-Gene Expression Correlations')




# 12. IDENTIFY HIGHLY CORRELATED TF-GENE PAIRS ---------------------------------

# Find TF-gene pairs with correlation > 0.8 (excluding self-correlations)
highCorPairs <- which(corMat > 0.8 & corMat != 1, arr.ind = TRUE)

# Number of highly correlated pairs
cat("Number of highly correlated TF-gene pairs:", nrow(highCorPairs), "\n")




# 13. ODDS RATIO ANALYSIS AT SINGLE-CELL LEVEL ---------------------------------

# Note: Pearson and Spearman correlations are problematic for scRNA-seq
# due to zero-inflation. Instead, use odds ratio to measure co-occurrence
# of TF and target gene expression in the same cell.
#
# Odds Ratio interpretation:
#   = 1: No association (expected co-occurrence)
#   > 1: Positive association (genes co-occur more than expected)
#   < 1: Negative association (genes co-occur less than expected)
#
# We use log odds ratio for visualization (log(1) = 0)

# Calculate odds ratios for highly correlated pairs
oddsRatios <- apply(highCorPairs,
                    1,
                    function(i) {
                      row <- rownames(corMat)[i[1]]
                      col <- colnames(corMat)[i[2]]
                      
                      # Count co-occurrence patterns
                      inBoth <- length(which(gbox[row, ] > 0 & gbox[col, ] > 0))
                      inNone <- length(which(gbox[row, ] == 0 & gbox[col, ] == 0))
                      inFirst <- length(which(gbox[row, ] > 0 & gbox[col, ] == 0))
                      inSecond <- length(which(gbox[row, ] == 0 & gbox[col, ] > 0))
                      
                      # Calculate odds ratio
                      oddsRatio <- (inBoth * inNone) / (inFirst * inSecond)
                      
                      c(inBoth, inFirst, inSecond, inNone, oddsRatio)
                    })

# Visualize distribution of log odds ratios
hist(log(oddsRatios[5, ]), 
     main = 'Log Odds Ratio Distribution - Flower',
     xlab = 'Log Odds Ratio',
     col = 'lightblue')




# 14. IDENTIFY STRONG TF-TARGET RELATIONSHIPS ----------------------------------

# Define threshold for "strong positive" relationships
# Using exp(1) ≈ 2.718 as the odds ratio threshold
thresh <- exp(1)

# Extract pairs with odds ratio > e (strong positive association)
strongPositiveIdx <- which(oddsRatios[5, ] > thresh)
strongPositive <- highCorPairs[strongPositiveIdx, ]

# Format for export (TF in column 1, target in column 2)
strongPositive[, 1] <- rownames(corMat)[strongPositive[, 1]]
strongPositive[, 2] <- colnames(corMat)[as.numeric(strongPositive[, 2])]
strongPositive <- cbind(strongPositive[, 2], 
                        strongPositive[, 1])  # TF comes before target
strongPositive <- cbind(strongPositive,
                        oddsRatios[5, strongPositiveIdx])

colnames(strongPositive) <- c("TF", "Target", "OddsRatio")

# Save strong positive regulatory pairs
write.table(strongPositive, 
            file = 'results/tables/Flower_doublePositives.txt', 
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

cat("Identified", nrow(strongPositive), "strong positive TF-target pairs\n")




# 15. IDENTIFY WEAK/NEGATIVE RELATIONSHIPS (SIMPSON'S PARADOX) -----------------

# Extract pairs with odds ratio < 1 (negative association at cell level)
# These are correlated at the cell type level but NOT at the single-cell level
# This represents Simpson's Paradox

thresh <- 1
simpsonIdx <- which(oddsRatios[5, ] < thresh) 
simpson <- highCorPairs[simpsonIdx, ]

# Format for export
simpson[, 1] <- rownames(corMat)[simpson[, 1]]
simpson[, 2] <- colnames(corMat)[as.numeric(simpson[, 2])]
simpson <- cbind(simpson[, 2], 
                 simpson[, 1]) 
simpson <- cbind(simpson,
                 oddsRatios[5, simpsonIdx])

colnames(simpson) <- c("TF", "Target", "OddsRatio")

# Save Simpson's paradox pairs
write.table(simpson,
            file = 'results/tables/Flower_SimpsonPairs.txt',
            sep = '\t',
            row.names = FALSE,
            quote = FALSE)

cat("Identified", nrow(simpson), "Simpson's paradox pairs\n")






# 16. NETWORK CONSTRUCTION WITH GENIE3 -----------------------------------------

# Additional packages required for network analysis
# BiocManager::install("GENIE3")
# install.packages('igraph')
# install.packages('network')
library(GENIE3)
library(igraph)
library(network)

# GENIE3 is a network inference algorithm based on random forests
# Parameters:
#   - Gene expression matrix (gbox)
#   - List of regulators (tfSubs - the G-box binding TFs)
#   - nTrees: Number of trees in random forest (reduced to 5 for speed)

# Construct gene regulatory network
net <- GENIE3(as.matrix(gbox), 
              regulators = tfSubs, 
              nTrees = 5)

# Save the network for future use
save(net,
     file = 'results/networks/Flower_Network_nTree_5.RData')

cat("Gene regulatory network constructed successfully\n")




# 17. CONVERT NETWORK TO EDGE LIST ---------------------------------------------

# Convert network matrix to adjacency list format
# Threshold of 0.05 keeps top-scoring regulatory relationships
ginieOutput <- convertToAdjacency(net, 0.05)

# Dimensions of the network (number of edges)
dim(ginieOutput)
cat("Network contains", nrow(ginieOutput), "regulatory edges\n")

# Preview first 10 edges
ginieOutput[1:10, ]




# 18. EXTRACT RANKED EDGE LIST -------------------------------------------------

# Get complete ranked list of regulatory edges
newNet <- GENIE3::getLinkList(net)

cat("Complete ranked network contains", nrow(newNet), "edges\n")




# 19. COMPARE NETWORKS: SINCERABOXCIS VS ARABOXCIS -----------------------------

# Compare the new single-cell network with the original bulk RNA-seq network

# Get unique genes present in the new network
genesInNet <- unique(c(newNet[, 1], newNet[, 2]))

# Filter AraBOXcis network to contain only genes present in new network
araboxcisFiltered <- araboxcis[which(araboxcis[, 1] %in% genesInNet & 
                                       araboxcis[, 2] %in% genesInNet), ]

# Extract top edges from new network (same number as filtered AraBOXcis)
newNetTopEdges <- newNet[1:nrow(araboxcisFiltered), ]

# Format edges for comparison
edgesNew <- paste(newNetTopEdges[, 1], newNetTopEdges[, 2], sep = '_')
edgesOld <- paste(araboxcisFiltered[, 1], araboxcisFiltered[, 2], sep = '_')

# Calculate overlaps
overlap <- length(which(edgesNew %in% edgesOld))
newOnly <- length(which(!(edgesNew %in% edgesOld)))
oldOnly <- length(which(!(edgesOld %in% edgesNew)))

cat("Network comparison:\n")
cat("  Edges in both networks:", overlap, "\n")
cat("  Edges only in SinceAraBOXcis:", newOnly, "\n")
cat("  Edges only in AraBOXcis:", oldOnly, "\n")




# 20. IDENTIFY IMPORTANT GENES: DEGREE ANALYSIS --------------------------------

# Degree = number of edges per transcription factor
# TFs with high degree regulate many downstream genes

# Count edges for each TF in both networks
tfsNew <- table(newNetTopEdges[, 1])
tfsOld <- table(araboxcisFiltered[, 1])[names(tfsNew)]

# Visualize TF degree distribution
hist(as.numeric(tfsNew), 
     main = 'TF Degree Distribution - SinceAraBOXcis', 
     xlab = 'Degree (Number of Target Genes)',
     col = 'lightblue')

# Compare TF degrees between networks
plot(as.numeric(tfsNew),
     as.numeric(tfsOld),
     xlab = 'Degree in SinceAraBOXcis',
     ylab = 'Degree in AraBOXcis',
     main = 'TF Degree Comparison Between Networks',
     pch = 20)

# Top 35 TFs by degree
cat("\nTop 35 TFs by degree:\n")
print(sort(tfsNew, decreasing = TRUE)[1:35])




# 21. NETWORK TOPOLOGY: CREATE IGRAPH OBJECT -----------------------------------

# Convert edge list to igraph network object for topology analysis
simple_network <- graph_from_edgelist(as.matrix(newNetTopEdges[, c(1, 2)]))

cat("Network created with", vcount(simple_network), "nodes and", 
    ecount(simple_network), "edges\n")




# 22. NODE BETWEENNESS CENTRALITY -----------------------------------------------

# Betweenness: Identifies genes that bridge different network regions
# High betweenness = gene lies on many shortest paths between other genes
# These genes are critical for information flow in the network

node_betweenness_all <- betweenness(simple_network)
node_betweenness <- node_betweenness_all[which(node_betweenness_all > 0)]

cat("\nTop 20 genes by betweenness centrality:\n")
print(sort(node_betweenness, decreasing = TRUE)[1:20])

# Visualize betweenness distribution
plot(sort(node_betweenness),
     main = 'Node Betweenness Distribution',
     xlab = 'Rank',
     ylab = 'Betweenness Centrality',
     pch = 20)




# 23. NODE ALPHA CENTRALITY -----------------------------------------------------

# Alpha centrality: Measures influence of a node and its neighbors
# High centrality = gene is well-connected and neighbors are well-connected
# Identifies influential regulatory hubs

node_centrality_all <- alpha_centrality(simple_network)
node_centrality <- node_centrality_all[which(node_centrality_all > 0)]

cat("\nTop 20 genes by alpha centrality:\n")
print(sort(node_centrality, decreasing = TRUE)[1:20])

# Visualize centrality distribution
plot(sort(node_centrality),
     main = 'Node Alpha Centrality Distribution',
     xlab = 'Rank',
     ylab = 'Alpha Centrality',
     pch = 20)




# 24. HUB SCORE ANALYSIS --------------------------------------------------------

# Hub score: Identifies genes that point to many important genes
# High hub score = gene regulates many well-connected targets
# Complementary to centrality measures

node_hub_all <- hub_score(simple_network)$vector
node_hub <- node_hub_all[which(node_hub_all > 0)]

cat("\nTop 20 genes by hub score:\n")
print(sort(node_hub, decreasing = TRUE)[1:20])

# Visualize hub score distribution
plot(sort(node_hub),
     main = 'Node Hub Score Distribution',
     xlab = 'Rank',
     ylab = 'Hub Score',
     pch = 20)




# 25. COMPARE NETWORK METRICS ---------------------------------------------------

# Compare betweenness vs alpha centrality
plot(node_betweenness_all, 
     node_centrality_all,
     xlab = 'Betweenness Centrality',
     ylab = 'Alpha Centrality',
     main = 'Betweenness vs Alpha Centrality',
     pch = 20)

# Compare hub score vs alpha centrality
plot(node_hub_all, 
     node_centrality_all,
     xlab = 'Hub Score',
     ylab = 'Alpha Centrality',
     main = 'Hub Score vs Alpha Centrality',
     pch = 20)




# 26. GO TERM ENRICHMENT ANALYSIS WITH PAFWAY ----------------------------------

# PAFway: Identifies Gene Ontology terms that are enriched upstream or 
# downstream of other GO terms in the network
# This reveals hierarchical relationships between biological processes

# Load functional annotation data
load('data/functionalData.RData')

# Run PAFway analysis
pafwayOut <- pafway(GOconcat, newNetTopEdges, unique(goOfInterest))
rownames(pafwayOut) <- colnames(pafwayOut) 

# Filter for rows/columns with at least one significant association (p < 0.05)
atLeastOneSigRow <- which(apply(pafwayOut,
                                1,
                                function(i) {
                                  length(which(i < 0.05))
                                }) > 0)

atLeastOneSigCol <- which(apply(pafwayOut,
                                2,
                                function(i) {
                                  length(which(i < 0.05))
                                }) > 0)

pafwayInterestingOnly <- pafwayOut[atLeastOneSigRow, atLeastOneSigCol]

cat("\nSignificant GO term associations found:", 
    length(atLeastOneSigRow), "terms\n")

# Visualize GO term hierarchy
# Column terms are upstream of row terms
# Smaller p-values = stronger associations
pheatmap(pafwayInterestingOnly,
         main = 'GO Term Network Hierarchy (p-values)')

# Log-transform for better visualization of highly significant associations
pheatmap(log10(pafwayInterestingOnly),
         main = 'GO Term Network Hierarchy (log10 p-values)')




# 27. SAVE NETWORK METRICS ------------------------------------------------------

# Save all centrality metrics for future analysis
save(node_betweenness, 
     node_centrality, 
     node_hub, 
     file = "results/tables/centrality_flower.RData")

save(node_betweenness_all, 
     node_centrality_all, 
     node_hub_all, 
     file = 'results/tables/centrality_all_flower.RData')

cat("\nAll network metrics saved successfully\n")




# ANALYSIS COMPLETE ------------------------------------------------------------
# Key outputs saved:
# - results/tables/Flower_avgExpressionByCluster.txt: Average gene expression per cell type
# - results/tables/Flower_doublePositives.txt: Strong TF-target relationships
# - results/tables/Flower_SimpsonPairs.txt: Paradoxical correlations (cell type vs single-cell)
# - results/networks/Flower_Network_nTree_5.RData: Inferred gene regulatory network
# - results/tables/centrality_flower.RData: Network centrality metrics (non-zero values)
# - results/tables/centrality_all_flower.RData: All network centrality metrics
#
# Network analysis complete:
# - GENIE3 network inference from single-cell data
# - Comparison with bulk RNA-seq AraBOXcis network
# - Identification of key regulatory hubs via multiple centrality metrics
# - GO term enrichment analysis revealing biological process hierarchy







