#' @title Align raw sample and group layers by shared identifiers.
#'
#' @description Filters raw sample and group layer data so that both contain the same identifiers in the same order before multi-omics networks are built.
#' Identifiers that occur only in the sample data or only in the group data are removed, because node IDs are generated positionally later in the workflow and must refer to the same biological identifiers in both networks.
#' Also removes identifiers containing "_" to match the LIONESS preprocessing step.
#'
#' @param sample_layers [list] Named list of sample-level layers.
#' @param group_layers [list] Named list of group-level layers.
#' @param connect_on [string] Name of the identifier column used for alignment, e.g. "gene_name".
#' @param min_overlap [float] Minimum fraction of identifiers that must be shared between sample and group data for each layer.
#'
#' @return [list] A named list with the aligned "sample_layers" and aligned "group_layers".
#'
#' @keywords internal
#' @noRd
align_sample_group_raw_layers <- function(sample_layers, group_layers, connect_on, min_overlap = 0.7) {

  # Iterate through every sample layer.
  for (layer_name in names(sample_layers)) {
    
    # Declare the raw data file for the sample and check if connect_on is a column in it..
    # Find out all identifiers (for example all gene_names, depending on what connect_on is).
    # Then search for any identifiers containing a '_' and drop them.
    sample_raw <- sample_layers[[layer_name]]$raw
    if (!connect_on %in% colnames(sample_raw)) {
      stop(paste0("Column '", connect_on, "' was not found in sample layer ", layer_name), call. = FALSE)
    }
    identifier <- as.character(sample_raw[[connect_on]])
    identifier_with_underscore <- grepl("_", identifier)
    if (any(identifier_with_underscore)) {
      message(paste0("Dropping unique identifier because their names contain '_' in sample layer '", layer_name, "': ", paste0(identifier[identifier_with_underscore], collapse = ", ")))
    }
    sample_raw <- sample_raw[!identifier_with_underscore, , drop = FALSE]
    
    # Now declare sample_ids and shared_ids, which will be compared later on, after shared_ids has been filled with all correct shared ids between the group and sample layers.
    sample_ids <- as.character(sample_raw[[connect_on]])
    shared_ids <- sample_ids
    
    # Iterate through every group.
    for (group_name in names(group_layers)) {
      
      # If the group has the same layer as the current sample_layer, this if statement is executed.
      if (layer_name %in% names(group_layers[[group_name]])) {
        
        # Declare the raw data file for the group and check if connect_on is a column in it.
        # Find out all identifiers (for example all gene_names, depending on what connect_on is).
        # Then search for any identifiers containing a '_' and drop them.
        g_raw <- group_layers[[group_name]][[layer_name]]$raw
        if (!connect_on %in% colnames(g_raw)) {
          stop(paste0("Column '", connect_on, "' was not found in group layer ", layer_name), call. = FALSE)
        }
        identifier <- as.character(g_raw[[connect_on]])
        identifier_with_underscore <- grepl("_", identifier)
        if (any(identifier_with_underscore)) {
          message(paste0("Dropping unique identifier because their names contain '_' in group '", group_name, "', layer '", layer_name, "': ", paste0(identifier[identifier_with_underscore], collapse = ", ")))
        }
        g_raw <- g_raw[!identifier_with_underscore, , drop = FALSE]
        group_layers[[group_name]][[layer_name]]$raw <- g_raw
        
        # shared_ids are now narrowed down to the shared identifiers between the current group layer and every identifier saved in shared_ids.
        shared_ids <- intersect(shared_ids, as.character(g_raw[[connect_on]]))
      }
    }
    
    # After shared_ids are filled with every group layer and sample layer overlap they are sorted.
    shared_ids <- sort(shared_ids)
    
    # The overlap between the identifiers saved in the sample layer and the real shared identifiers between sample and group layers is calculated and evaluated, if too low the code stops and warns.
    # If the overlap check passes the sample layer is updated with the sorted shared_ids as identifier column, because of the nature of match, duplicates are also dropped here.
    sample_overlap <- length(shared_ids) / length(sample_ids)
    if (sample_overlap < min_overlap) {
      stop(paste0("Too little identifier overlap in layer '", layer_name, "' for sample data. Sample overlap: ", round(sample_overlap * 100, 1), "%."), call. = FALSE)
    }
    sample_layers[[layer_name]]$raw <- sample_raw[match(shared_ids, sample_ids), , drop = FALSE]
    
    # The group layers are also updated with the sorted shared_ids as identifier column, because of the nature of match, duplicates are also dropped here.
    # For this the code iterates through every group again.
    for (group_name in names(group_layers)) {
      
      # Any layer in the group that is the same as the sample layer we currently look at is catched here.
      if (layer_name %in% names(group_layers[[group_name]])) {
        
        # Declare the raw data file for the group.
        # Find out all identifiers (for example all gene_names, depending on what connect_on is).
        g_raw <- group_layers[[group_name]][[layer_name]]$raw
        identifier <- as.character(g_raw[[connect_on]])
        
        # The overlap between the identifiers saved in the group layer and the real shared identifiers between sample and group layers is calculated and evaluated, if too low the code stops and warns.
        # If the overlap check passes the group layer is updated with the sorted shared_ids as identifier column, because of the nature of match, duplicates are also dropped here.
        group_overlap <- length(shared_ids) / length(identifier)
        if (group_overlap < min_overlap) {
          stop(paste0("Too little identifier overlap in layer '", layer_name, "' for group '", group_name, "'. Group overlap: ", round(group_overlap * 100, 1), "%."), call. = FALSE)
        }
        group_layers[[group_name]][[layer_name]]$raw <- g_raw[match(shared_ids, identifier), , drop = FALSE]
      }
    }
  }
  
  return(list(
    sample_layers = sample_layers,
    group_layers  = group_layers
  ))
}

#' @title Prepare raw molecular data for LIONESS.
#'
#' @description Formats a raw molecular data table into the matrix format required by LIONESS.
#' The non numerical metadata columns are removed, the identifier column is used as row names and identifiers containing "_" are dropped because LIONESS interaction names are later split by this character.
#' Also checks if identifier column is unique and if not drops duplicates and prints the duplications.
#'
#' @param raw_data [data.frame] Raw data. Rows represent unique identifiers such as genes or proteins, columns contain metadata and sample measurements.
#' @param identifier_col [string] Name of the column containing the unique identifier, e.g. "gene_name".
#' @param meta_cols [character vector] Names of metadata columns that have non numerical data and do not contain sample specific measurements.
#'
#' @return A matrix with identifiers as row names and samples as columns.
#'
#' @keywords internal
#' @noRd
prepare_lioness_expr <- function(raw_data, identifier_col, meta_cols) {

  # Check if the identifier column really exists in the raw data.
  # Find out all identifiers (for example all gene_names, depending on what connect_on is).
  # Then search for any identifiers containing a '_' and drop them.
  if (!identifier_col %in% colnames(raw_data)) {
    stop(paste0("Column '", identifier_col, "' was not found in the raw data of the samples."), call. = FALSE)
  }
  identifier <- as.character(raw_data[[identifier_col]])
  identifier_with_underscore <- grepl("_", identifier)
  if (any(identifier_with_underscore)) {
    message(paste0("Dropping unique identifier because their names contain '_': ", paste0(identifier[identifier_with_underscore], collapse = ", ")))
  }
  raw_data <- raw_data[!identifier_with_underscore, , drop = FALSE]
  
  # Now search for any duplicates in the new identifiers and drop them if they exist.
  identifier <- as.character(raw_data[[identifier_col]])
  if (length(unique(identifier)) != length(identifier)) {
    duplicate_ids <- unique(identifier[duplicated(identifier)])
    message(paste0("Dropping identifiers ", paste0(duplicate_ids, collapse = ", "), " because they occur multiple times"))
  }
  raw_data <- raw_data[match(unique(identifier), identifier), , drop = FALSE]
  
  # Then the meta_cols are deleted, the new variable expr is changed to a matrix and gets the rownames of the identifier column without '_' or duplicates.
  expr  <- raw_data[, setdiff(colnames(raw_data), meta_cols)]
  expr <- as.matrix(expr)
  rownames(expr) <- as.character(raw_data[[identifier_col]])
  return(expr)
}

#' @title Generate single-sample networks with LIONESS.
#'
#' @description Runs LIONESS on one layer and converts the resulting assay matrix into a list of sample-specific edge tables.
#' Each list element contains the interaction name (for example GeneA_GeneB) in one column and the corresponding LIONESS edge weight for one sample in the other column.
#' The aggregate network is a Pearson correlation computed with WGCNA::cor, the same way DrDimont computes its correlation matrices, so missing values are handled as in DrDimont.
#'
#' @param expr [matrix] Matrix with identifiers in rows and samples in columns.
#' @param handling_missing_data ["all.obs"|"pairwise.complete.obs"] Specifying the handling of missing data during correlation computation, like DrDimont's handling_missing_data setting. (default: all.obs)
#'
#' @return A named list of data.frames. Each list element corresponds to one sample and contains the columns "interaction" and "weight".
#'
#' @keywords internal
#' @noRd
run_lioness <- function(expr, handling_missing_data = "all.obs") {

  # The aggregate network function is LIONESS' default (Pearson correlation), but computed with WGCNA::cor and handling_missing_data like in DrDimont.
  network_function <- function(x) {
    WGCNA::cor(t(x), method = "pearson", use = handling_missing_data)
  }

  # Lioness is run and results are written into lioness_out, then the results are extracted to lioness_mat.
  lioness_out <- lionessR::lioness(expr, f = network_function)
  lioness_mat <- SummarizedExperiment::assay(lioness_out, "lioness")
  
  # A list for all the single sample molecular networks is created and filled by looping through every lioness result and writing it into the new list.
  sample_list <- list()
  for (sample in colnames(lioness_mat)) {
    df <- data.frame(
      interaction = rownames(lioness_mat),
      weight      = lioness_mat[, sample],
      row.names   = NULL
    )
    sample_list[[sample]] <- df
  }
  return(sample_list)
}
