#' @title Pre-compute interaction scores for a single network.
#'
#' @description Runs the drug-target identification and interaction score graph generation for a single DrDimont-compatible multi-omics network in isolation.
#' The result can be cached and reused across multiple comparisons that involve this network.
#'
#' @param network [list] DrDimont-compatible multi-omics network object.
#' @param drug_target_interactions [list] DrDimont drug-target interaction object.
#' @param settings [list] DrDimont settings object.
#'
#' @return A named list with:
#' \describe{
#'   \item{interaction_score_graph}{The interaction score graph for this network (from "generate_interaction_score_graphs"$groupA).}
#'   \item{drug_targets}{The drug target mapping for this network (from "determine_drug_targets"$targets).}
#' }
#'
#' @keywords internal
#' @noRd
compute_single_network_scores <- function(network, drug_target_interactions, settings) {
  # The describe \item part of the comments has been taken by suggestion of ChatGPT (LLM2)

  # Save the multi-omics graph and the annotations in two different variables.
  single_graphs <- list(groupA = network$graphs$groupA)
  single_annotations <- list(both = network$annotations$both)
  
  # Calculate the drug target edges inside the multi-omics graph with the help of the DrDimont function determine_drug_targets and save them in drug_target_edges.
  drug_target_edges <- DrDimont::determine_drug_targets(
    graphs = single_graphs,
    annotations = single_annotations,
    drug_target_interactions = drug_target_interactions,
    settings = settings
  )
  
  # Compute the interaction score graph with the help of the DrDimont function generate_interaction_score_graphs and store it in interaction_score_graphs.
  interaction_score_graphs <- DrDimont::generate_interaction_score_graphs(
    graphs = single_graphs,
    drug_target_edgelists = drug_target_edges[["edgelists"]],
    settings = settings
  )
  
  return(list(
    interaction_score_graph = interaction_score_graphs$groupA,
    drug_targets = drug_target_edges[["targets"]]
  ))
}

#' @title Compare two DrDimont-compatible multi-omics networks for their drug interations.
#'
#' @description Computes differential drug response scores between two DrDimont-compatible multi-omics network objects.
#'
#' @param drug_target_interactions [list] DrDimont drug-target interaction object.
#' @param settings [list] DrDimont settings object.
#' @param precomputed_a [list] Pre-computed scores from "compute_single_network_scores".
#' @param precomputed_b [list] Pre-computed scores from "compute_single_network_scores".
#'
#' @return A data.frame containing drug names and their differential drug response scores.
#'
#' @keywords internal
#' @noRd
compare_drug_interactions <- function(drug_target_interactions, settings, precomputed_a, precomputed_b) {

  # Store the interaction score graphs to be compared to each other in the variable interaction_score_graphs so that DrDimont recognizes them later on.
  interaction_score_graphs <- list(
    groupA = precomputed_a$interaction_score_graph,
    groupB = precomputed_b$interaction_score_graph
  )
  
  # Create a variable named drug_targets to save the drug targets from the interactions score graphs.
  # Because they should be the same, here precomputed_a is used to look into them.
  drug_targets <- precomputed_a$drug_targets
  
  # Compute a differential graph between the two interaction score graphs with the help of the DrDimont function generate_differential_score_graph.
  differential_graph <- DrDimont::generate_differential_score_graph(
    interaction_score_graphs = interaction_score_graphs,
    settings = settings
  )
  
  # Calculate the drug response scores for the differentail graph with the help of the DrDimont function compute_drug_response_scores.
  drug_response_scores <- DrDimont::compute_drug_response_scores(
    differential_graph = differential_graph,
    drug_targets = drug_targets,
    settings = settings
  )

  return(drug_response_scores)
}

#' @title Validate identifier order between sample and group annotations.
#'
#' @description Checks whether the identifiers in the "connect_on" column occur in the same order between a sample annotation table and a group annotation table for one layer.
#' This is required because graph vertices are named by node_id, and the same node_id must refer to the same biological identifier in both graphs.
#'
#' @param sample_ann [data.frame] Annotation table of a sample-level multi-omics network.
#' @param group_ann [data.frame] Annotation table of a group-level multi-omics network.
#' @param connect_on [string] Name of the identifier column to compare, e.g. "gene_name".
#' @param layer_name [string] Name of the layer for which identifiers should be compared.
#'
#' @return No return value.
#'
#' @keywords internal
#' @noRd
validate_connect_on_alignment <- function(sample_ann, group_ann, connect_on, layer_name) {

  # Save the annotations for the specific layer for the sample and group annotations in two variables.
  sample_layer_ann <- sample_ann[sample_ann$layer == layer_name, , drop = FALSE]
  group_layer_ann  <- group_ann[group_ann$layer == layer_name, , drop = FALSE]
  
  # Save the identifiers into two new variables for samples and groups.
  sample_vals <- as.character(sample_layer_ann[[connect_on]])
  group_vals  <- as.character(group_layer_ann[[connect_on]])
  
  # Check if the identifiers are equal length for the sample and group data, if not, stop the function.
  if (length(sample_vals) != length(group_vals)) {
    stop(paste0("Different number of identifiers in layer '", layer_name, "'.\n", "Sample has ", length(sample_vals), " identifiers.\n", "Group has ", length(group_vals), " identifiers."), call. = FALSE)
  }
  
  # Check if the values of sample and group are the same in every position, if not, stop the function.
  if (!identical(sample_vals, group_vals)) {
    mismatch_id <- which(sample_vals != group_vals)
    first_mismatches <- utils::head(mismatch_id, 5)
    mismatch_text <- paste(paste0("position ", first_mismatches, ": sample = ", sample_vals[first_mismatches], ", group = ", group_vals[first_mismatches]), collapse = "\n")
    stop(paste0("Identifier order differs in layer '", layer_name, "'.\n", mismatch_text), call. = FALSE)
  }
  
  message("[OK] Layer '", layer_name, "': connect_on values match in the same order (", length(sample_vals), " entries).")
}