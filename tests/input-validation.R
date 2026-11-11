# Checks that wrong input to create_lioness_multiomics_samples() is caught by
# check_lioness_input() (R/checks.R) with clear messages, all reported together,
# before any network is built.
# R CMD check runs this script; it fails as soon as a stopifnot() condition is not met.

library(SingleSampleDrDimont)
check_lioness_input <- SingleSampleDrDimont:::check_lioness_input

mrna <- data.frame(gene_name = c("A", "B", "C"), s1 = c(1, 2, 3), s2 = c(2, 1, 4))
valid_layers <- list(mrna = list(raw = mrna, meta_cols = "gene_name"))

# Valid input gives no errors.
stopifnot(length(check_lioness_input(valid_layers, NULL, "gene_name")) == 0)

# A missing connect_on gives one error that names it.
errors <- check_lioness_input(valid_layers, NULL, NULL)
stopifnot(length(errors) == 1, grepl("connect_on", errors))

# Two mistakes in one layer are both reported: a meta column that is not in raw,
# and text in a measurement column.
broken <- data.frame(gene_name = c("A", "B", "C"), s1 = c("x", "y", "z"))
errors <- check_lioness_input(
  list(mrna = list(raw = broken, meta_cols = c("gene_name", "ref_seq"))), NULL, "gene_name"
)
stopifnot(
  length(errors) == 2,
  any(grepl("ref_seq", errors)),
  any(grepl("not numeric", errors))
)

# The workflow stops on these errors before running LIONESS.
result <- tryCatch(
  suppressMessages(create_lioness_multiomics_samples(valid_layers, connect_on = NULL)),
  error = function(e) e
)
stopifnot(inherits(result, "error"))
