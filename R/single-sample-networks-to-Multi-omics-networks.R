#' @title Convert LIONESS edge table to adjacency matrix.
#'
#' @description Converts a LIONESS edge table into an adjacency matrix.
#' Edge weights outside the interval [-1, 1] are either set to zero or transformed with tanh(), because DrDimont functions expect edge weights in this range.
#' Undefined edge weights (NaN) are set to zero, i.e. no edge.
#'
#' @param sample [data.frame] LIONESS edge table containing the columns "interaction" (in the form of "identifier1_identifier2") and "weight".
#' @param clip_LIONESS_weights [logical] If "TRUE", weights outside [-1, 1] are transformed with tanh().
#' If "FALSE", weights outside [-1, 1] are set to zero.
#'
#' @return An adjacency matrix with identifiers as row and column names.
#'
#' @keywords internal
#' @noRd
edge_to_adjacency <- function(sample, clip_LIONESS_weights) {

  # The interactions from the Lioness results are split at the '_' (for example "GeneA_GeneB" gets split into c("GeneA", "GeneB")).
  parts <- strsplit(sample$interaction, split = "_")
  
  # A data frame named x is created which contains the necessary values for building the adjacency matrix.
  x <- data.frame(
    interaction = sample$interaction,
    weight = sample$weight,
    identifier1 = vapply(parts, `[`, character(1), 1),
    identifier2 = vapply(parts, `[`, character(1), 2),
    stringsAsFactors = FALSE
  )

  # Weights are undefined (NaN) when an identifier has no variance in a leave-one-out network of LIONESS, for example when it is measured at the same value in every other sample.
  # These weights are set to 0 (no edge), because there is no correlation to describe, and the user is informed.
  undefined_weight <- is.na(x$weight)
  if (any(undefined_weight)) {
    message(paste0("Setting ", sum(undefined_weight), " undefined LIONESS edge weights (NaN) to 0. This happens when an identifier has no variance in a leave-one-out network."))
    x$weight[undefined_weight] <- 0
  }

  # Create an identifier variable which contains every identifier once.
  # Then these identifiers are used to build a matrix containing all identifiers horizontally and vertically, with all values between them being 0.
  identifier <- unique(c(x$identifier1, x$identifier2))
  mat <- matrix(0, nrow = length(identifier), ncol = length(identifier), dimnames = list(identifier, identifier))
  
  # The matrix is filled with the corresponding weights for the identifiers.
  # If the weights are outside the range of [-1, 1], they are either clipped between [-1, 1] via tanh(), or set to 0, depending on what clip_LIONESS_weights is.
  for (i in seq_len(nrow(x))) {
    i1 <- x$identifier1[i]
    i2 <- x$identifier2[i]
    w <- x$weight[i]
    if (clip_LIONESS_weights) {
      if (w < -1 || w > 1) {
        w <- tanh(w)
      }
    } else {
      if (w < -1 || w > 1) {
        w <- 0
      }
    }
    mat[i1, i2] <- w
    mat[i2, i1] <- w
  }
  return(mat)
}

#' @title Create DrDimont-compatible annotation table
#'
#' @description Creates an annotation table that maps node IDs to identifiers.
#'
#' @param node_names [character vector] Row names of the adjacency matrix (the identifiers).
#' @param layer_name [string] Name of the molecular layer.
#' @param identifier_col [string] Name of the identifier column to create in the annotation table, e.g. "gene_name".
#'
#' @return A data.frame with the columns "node_id", "identifier_col", and "layer".
#'
#' @keywords internal
#' @noRd
make_annotation <- function(node_names, layer_name, identifier_col) {

  # Create ann variable, which stores node_id (for example for some sample from the mrna layer and its first row node_id would be: mrna_1) and the layer name.
  # The identifier column gets added too, but not inside data.frame because then identifier_col would be read as literal "identifier_col" and not as what it is as a variable, like "gene_name" for example.
  ann <- data.frame(
    node_id = paste0(layer_name, "_", seq_along(node_names)),
    layer = layer_name,
    stringsAsFactors = FALSE
  )
  ann[[identifier_col]] <- node_names
  return(ann)
}

#' @title Convert adjacency matrix to igraph object.
#'
#' @description Converts a adjacency matrix into an undirected igraph object and replaces the original vertex names with DrDimont-compatible node IDs (e.g. mrna_5, protein_4, depending on layer and number of row).
#' The resulting graph has all self-connections set to zero.
#'
#' @param adj_mat [matrix] Weighted adjacency matrix.
#' @param annotation [data.frame] Annotation table containing at least "node_id" and the identifier column specified by "identifier_col".
#' @param identifier_col [string] Name of the identifier column used to match matrix row names to annotation rows, e.g. "gene_name"..
#'
#' @return An undirected weighted igraph object with vertex names set to "node_id".
#'
#' @keywords internal
#' @noRd
adj_to_igraph <- function(adj_mat, annotation, identifier_col) {

  # Create igraph g.
  # diag = FALSE ensures there are no self connections.
  g <- igraph::graph_from_adjacency_matrix(
    adjmatrix = adj_mat,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )

  # The variable id looks at the names of the graph nodes and matches to the annotations, so it contains the positions where in the annotation table the graph nodes are.
  # It looks like c(2, 3, 1) if the graph nodes 1, 2, 3 are at position 2, position 3 and position 1 in the annotation table).
  # Then the graph nodes are named like the node_id is named for that position, so at the start the name of the graph nodes was for example GeneA, but now its mrna_1.
  # This was done to ensure compatibility with DrDimont in the later comparison because the graphs from DrDimont have these naming conventions.
  id <- match(igraph::V(g)$name, annotation[[identifier_col]])
  igraph::V(g)$name <- as.character(annotation$node_id[id])
  return(g)
}

#' @title Build a DrDimont-compatible multi-omics network for one sample.
#'
#' @description Builds one combined multi-omics graph from a list of single-sample layer-specific adjacency matrices.
#' For each layer, annotations are created, adjacency matrices are converted to igraph objects, and inter-layer connections are generated using the identifier column specified by "identifier_col".
#'
#' @param sample_adjs [list] Named list of adjacency matrices for one sample.
#' Each list element corresponds to one layer.
#' @param settings [list] DrDimont settings object.
#' @param identifier_col [string] Name of the identifier column used for inter-layer connections, e.g. "gene_name".
#'
#' @return A DrDimont-compatible multi-omics network object containing the sample graph and its node annotations.
#'
#' @keywords internal
#' @noRd
build_sample_multilayer <- function(sample_adjs, settings, identifier_col) {

  # Create a variable containing all the layer names of the sample.
  layer_names <- names(sample_adjs)
  
  # Sort all rows and columns for the matrices.
  # This is to ensure later that they are in the same order as the group networks.
  sample_adjs <- lapply(sample_adjs, function(mat) {
    ord <- sort(rownames(mat))
    mat[ord, ord, drop = FALSE]
  })
  
  # Create annotations variable to hold layer specific annotations.
  # Then populate the variable with these annotations for each layer.
  annotations_by_layer <- list()
  for (layer_name in layer_names) {
    ann <- make_annotation(
      node_names = rownames(sample_adjs[[layer_name]]),
      layer_name = layer_name,
      identifier_col = identifier_col
    )
    annotations_by_layer[[layer_name]] <- ann
  }
  
  # Create a variable to hold all graphs for each layer.
  # Then fill it with graphs for each layer.
  graphs_by_layer <- list()
  for (layer_name in layer_names) {
    graphs_by_layer[[layer_name]] <- adj_to_igraph(
      adj_mat = sample_adjs[[layer_name]],
      annotation = annotations_by_layer[[layer_name]],
      identifier_col = identifier_col
    )
  }
  
  # Create a graphs and annotations object for DrDimont function generate_combined_graphs.
  # Because only single sample networks exist, groupB is left empty, for only one group exists.
  # annotations$both needs to contain unique identifiers across the whole data, so the layer annotations get combined here to fit the schema.
  graphs <- list(
    groupA = graphs_by_layer,
    groupB = NULL
  )
  annotations <- list(
    groupA = annotations_by_layer,
    groupB = NULL,
    both = do.call(rbind, annotations_by_layer)
  )
  
  # Create a variable containing all possible combinations for layers, so multi-omics networks get created not only for a single possible layer combination for the sample, but all of them.
  inter_layer_pairs <- utils::combn(layer_names, 2, simplify = FALSE)
  
  # Creates the list inter_layer_connections that specifies the connection between two layers, every pair from inter_layer_pairs is used for this.
  # The connections between the two layers are calculated via the DrDimont function make_connection.
  inter_layer_connections <- lapply(inter_layer_pairs, function(pair) {
    DrDimont::make_connection(
      from = pair[[1]],
      to = pair[[2]],
      connect_on = identifier_col,
      weight = 1,
      group = "A"
    )
  })
  
  # Combines the layers of the single sample together into one graph with the hellp of the DrDimont function generate_combined_graphs and saves them in the variable combined.
  combined <- DrDimont::generate_combined_graphs(
    graphs = graphs,
    annotations = annotations,
    inter_layer_connections = inter_layer_connections,
    settings = settings
  )
  
  return(combined)
}

#' @title Save intermediate workflow results to directory.
#'
#' @description Saves an R object as an ".rds" file.
#'
#' @param object [object] The R object to save.
#' @param filename [string] Name of the ".rds" file to create inside "intermediate_dir".
#' @param dir [string] Directory path for saved intermediate results.
#'
#' @return No return value.
#'
#' @keywords internal
#' @noRd
result_saver <- function(object, filename, dir) {

  # Saves an R object with a specified filename to the specified path.
  path <- file.path(dir, filename)
  saveRDS(object, path)
  message("Saved to: ", path)
}