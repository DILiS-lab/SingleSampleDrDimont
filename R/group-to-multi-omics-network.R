#' @title Build a DrDimont-compatible multi-omics network for one group.
#'
#' @description Builds a group-level multi-omics network from raw data using the standard DrDimont graph generation workflow.
#' Each raw layer is converted with "make_layer", correlation matrices are computed and individual layer graphs are generated.
#' Then, the resulting graphs are combined using inter-layer connections based on "connect_on".
#' If the reduction with "pickHardThreshold" fails, the layer graphs are generated again with DrDimont's "p_value" reduction.
#'
#' @param group_layers [list] Named list of layers for one group.
#' Each element must contain "raw", a data.frame with measurements and metadata, and "meta_cols", a character vector naming metadata columns with non numerical data.
#' @param settings [list] DrDimont settings object.
#' @param connect_on [string] Name of the identifier column used to connect nodes between layers, e.g. "gene_name".
#'
#' @return [list] A DrDimont-compatible multi-omics network object containing the combined group graph and its node annotations.
#'
#' @keywords internal
#' @noRd
build_group_multilayer_drdimont <- function(group_layers, settings, connect_on) {

  # Sort the raw data for the layers according to their identifier column.
  group_layers <- lapply(group_layers, function(layer) {
    layer$raw <- layer$raw[order(layer$raw[[connect_on]]), , drop = FALSE]
    return(layer)
  })
  
  # Create a new variable called drdimont_layers and fill it with layers computed with DrDimont for every layer passed in by the user.
  drdimont_layers <- lapply(names(group_layers), function(layer_name) {
    
    # Create some variables, containing the current layer name, the raw data for this layer and meta columns for this layer, respectively.
    layer <- group_layers[[layer_name]]
    raw_data <- layer$raw
    meta_cols <- layer$meta_cols
    
    # Save all the column containg numerical data inside the variable data_cols.
    # Then use the data to build a layer with the help of the DrDimont function make_layer.
    data_cols <- setdiff(colnames(raw_data), meta_cols)
    DrDimont::make_layer(
      name = layer_name,
      data_groupA = raw_data[, data_cols, drop = FALSE],
      data_groupB = NULL,
      identifiers_groupA = raw_data[, meta_cols, drop = FALSE],
      identifiers_groupB = NULL
    )
  })
  
  # Name the layers in DrDimont after the names they had originally from the user.
  # Then save these names in a variable and create a variable containing all possible combinations for layers, so multi-omics networks get created not only for a single possible layer combination for the group, but all of them.
  names(drdimont_layers) <- names(group_layers)
  layer_names <- names(drdimont_layers)
  inter_layer_pairs <- utils::combn(layer_names, 2, simplify = FALSE)
  
  # Creates the list inter_layer_connections that specifies the connection between two layers, every pair from inter_layer_pairs is used for this.
  # The connections between the two layers are calculated via the DrDimont function make_connection.
  inter_layer_connections <- lapply(inter_layer_pairs, function(pair) {
    DrDimont::make_connection(
      from = pair[[1]],
      to = pair[[2]],
      connect_on = connect_on,
      weight = 1,
      group = "A"
    )
  })
  
  # Generate correlation matrices necessary for computing the individual graphs for the group layers.
  correlation_matrices <- DrDimont::compute_correlation_matrices(
    layers = drdimont_layers,
    settings = settings
  )
  
  # For cuts above the strongest correlation of a layer no node has any edge, and WGCNA's pickHardThreshold then stops in browser() instead of throwing an error.
  # So if the cut_vector of any pickHardThreshold layer goes above its strongest correlation, the group network uses DrDimont's p_value reduction instead.
  pickHardThreshold_layers <- Filter(function(layer_name) {
    get_layer_setting(layer_name, "", settings, "reduction_method") == "pickHardThreshold"
  }, names(drdimont_layers))
  cut_above_max_correlation <- Filter(function(layer_name) {
    adjacency_matrix <- correlation_matrices$correlation_matrices$groupA[[layer_name]]
    max(get_layer_setting(layer_name, "groupA", settings, "cut_vector")) > max(abs(adjacency_matrix[upper.tri(adjacency_matrix)]), na.rm = TRUE)
  }, pickHardThreshold_layers)
  if (length(cut_above_max_correlation) > 0) {
    message("cut_vector goes above the strongest correlation in layer(s) ", paste0(cut_above_max_correlation, collapse = ", "), ", so pickHardThreshold cannot reduce the group network. Falling back to p_value reduction.")
    settings$reduction_method <- "p_value"
  }

  # Generate individual_graphs necessary for computing a combined_group_graph for the group layers.
  # If the reduction with pickHardThreshold fails otherwise (for example when no cut reaches the R^2 cutoff), the graphs are generated again with DrDimont's p_value reduction.
  individual_graphs <- tryCatch(
    DrDimont::generate_individual_graphs(
      correlation_matrices = correlation_matrices,
      layers = drdimont_layers,
      settings = settings
    ),
    error = function(e) {
      if (!"pickHardThreshold" %in% unlist(settings$reduction_method)) {
        stop(e)
      }
      message("pickHardThreshold reduction of the group network failed (", conditionMessage(e), "). Falling back to p_value reduction.")
      settings$reduction_method <- "p_value"
      DrDimont::generate_individual_graphs(
        correlation_matrices = correlation_matrices,
        layers = drdimont_layers,
        settings = settings
      )
    }
  )
  
  # Combines the layer graphs of the group together into one graph with the help of the DrDimont function generate_combined_graphs and saves them in the variable combined_group_graph.
  combined_group_graph <- DrDimont::generate_combined_graphs(
    graphs = individual_graphs$graphs,
    annotations = individual_graphs$annotations,
    inter_layer_connections = inter_layer_connections,
    settings = settings
  )
  
  return(combined_group_graph)
}