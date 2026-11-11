# Checks how LIONESS edge weights become the adjacency matrices DrDimont works with
# (edge_to_adjacency() in R/single-sample-networks-to-Multi-omics-networks.R):
# undefined weights (NaN) become 0, and weights outside [-1, 1] become tanh(weight)
# or 0, depending on clip_LIONESS_weights.
# R CMD check runs this script; it fails as soon as a stopifnot() condition is not met.

library(SingleSampleDrDimont)
edge_to_adjacency <- SingleSampleDrDimont:::edge_to_adjacency

edges <- data.frame(
  interaction = c("A_B", "A_C", "B_C"),
  weight = c(0.5, 2, NaN)
)

clipped <- suppressMessages(edge_to_adjacency(edges, clip_LIONESS_weights = TRUE))
stopifnot(
  isSymmetric(clipped),
  all(diag(clipped) == 0),
  clipped["A", "B"] == 0.5,
  clipped["A", "C"] == tanh(2),
  clipped["B", "C"] == 0
)

zeroed <- suppressMessages(edge_to_adjacency(edges, clip_LIONESS_weights = FALSE))
stopifnot(
  zeroed["A", "B"] == 0.5,
  zeroed["A", "C"] == 0,
  zeroed["B", "C"] == 0
)
