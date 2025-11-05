

# 1. Export the flower network edge list
# This uses your GENIE3 network 
write.csv(newNetTopEdges, file = "~/Documents/Bioinformatics/Araboxcis/data/exported_data/Flower_network.csv", row.names = FALSE)


# 2. Export the flower GO enrichment results (the pafway output)
# First create a long-format table from your pafwayInterestingOnly matrix
go_results <- data.frame()
for (row_term in rownames(pafwayInterestingOnly)) {
  for (col_term in colnames(pafwayInterestingOnly)) {
    go_results <- rbind(go_results, data.frame(
      DownstreamTerm = row_term,
      UpstreamTerm = col_term,
      pValue = pafwayInterestingOnly[row_term, col_term]
    ))
  }
}

# Export the GO results
write.csv(go_results, file = "~/Documents/Bioinformatics/Araboxcis/data/exported_data/Flower_GO_results.csv", row.names = FALSE)


# 3. Export the GO terms that were found to be significant 
sig_go_terms <- unique(c(rownames(pafwayInterestingOnly), colnames(pafwayInterestingOnly)))
write.csv(data.frame(GOTerm = sig_go_terms), file = "~/Documents/Bioinformatics/Araboxcis/data/exported_data/Flower_sig_GO_terms.csv", row.names = FALSE)

# Optional: Export the full pafwayOut matrix before filtering for significant terms
# This contains all GO terms analysed, not just the significant ones
write.csv(pafwayOut, file = "~/Documents/Bioinformatics/Araboxcis/data/exported_data/Flower_full_pafway_matrix.csv")
