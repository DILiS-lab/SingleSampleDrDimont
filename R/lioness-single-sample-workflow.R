#' @title Create LIONESS-based single-sample multi-omics networks and compare them if needed.
#'
#' @description Runs a single-sample multi-omics workflow based on LIONESS and DrDimont.
#' The raw data files for the layers are first formatted for LIONESS.
#' Single-sample networks are generated for each layer and are converted to adjacency matrices.
#' Then they get reduced with DrDimont's "pickHardThreshold" network reduction and combined into DrDimont-compatible multi-omics networks.
#' Optionally, two samples can be compared or one sample can be compared to a group-level DrDimont network.
#'
#' @param sample_layers [list] Named list of layers used to create single-sample networks.
#' Each element must contain "raw", a data.frame containing all the measured information.
#' Each element must also contain "meta_cols", a character vector naming the non numerical metadata columns.
#' All other columns must be numeric.
#' @param group_layers [list|null] Optional named list of group-level layers. Required when "compare$type" is "sample_vs_group".
#' Each element must contain "raw" and "meta_cols" like sample_layers.
#' @param drug_target_interactions [list|null] Optional DrDimont drug-target interaction object.
#' Required when "compare" is not "NULL" and the goal is to compute not only multi-omics networks but also drug interactions.
#' @param compare [list|null] Optional comparison specification.
#' For sample-vs-sample comparison, use "list(type = "sample_vs_sample", samples = list(c("sample_1", "sample_2")))", with "sample_1" and "sample_2" being the names of the columns for the samples to be compared.
#' For sample-vs-group comparison, use "list(type = "sample_vs_group", samples = list(c("sample_1", "group1))",  with "sample_1" being the name of the column for the sample to be compared.
#' @param r_squared_cutoff_samples [float] Minimum scale-free topology fitting index.
#' Used for reducing LIONESS-derived sample networks with the DrDimont method "network_reduction_by_pickHardThreshold".
#' @param r_squared_cutoff_groups [float] Minimum scale-free topology fitting index.
#' Used for reducing group-level DrDimont networks when "compare$type" is "sample_vs_group" with the DrDimont method "network_reduction_by_pickHardThreshold".
#' @param settings [list] DrDimont settings object.
#' Its "handling_missing_data" setting is also used for the LIONESS sample networks (see \code{\link[DrDimont]{drdimont_settings}}).
#' @param connect_on [string] Name of the identifier column used for annotation, inter-layer connections, and drug-target matching, e.g. "gene_name".
#' Required. It has to be listed in "meta_cols" of every layer.
#' @param clip_LIONESS_weights [logical] If "TRUE", weights outside [-1, 1] are transformed with tanh().
#' If "FALSE", weights outside [-1, 1] are set to zero.
#' @param min_layer [int] Defines the minimum amount of layers in which a sample has to occur to build a multi omics network with it.
#' @param min_overlap [float] Defines the minimum percentage of how much overlap should be between the passed genes in the raw data files for the sample and the group.
#' @param save_intermediate [logical] Whether to save intermediate multi-omics network objects as ".rds" files while the workflow is running.
#' @param intermediate_dir [string|null] Directory path for saved intermediate results. Required if "save_intermediate" is "TRUE"; created automatically if it does not exist.
#' @param save_results [logical] Whether to save multi-omics network objects as ".rds" files while the workflow is running.
#' @param results_dir [string|null] Directory path for saved results. Required if "save_results" is "TRUE"; created automatically if it does not exist.
#' @param group_reduction_method [string] Defines the method for group network reduction, the only valid inputs are "pickHardThreshold" and "p_value". Default is "pickHardThreshold".
#' If "pickHardThreshold" fails for a group network, that group network is reduced with "p_value" instead.
#'
#' @return A named list with the elements multiomics_per_sample, interaction_score_cache, group_multiomics and comparison_result
#'
#' @examples
#' data(mrna_data, package = "DrDimont")
#' data(protein_data, package = "DrDimont")
#'
#' # Use a subset of samples common to both omics layers, for a fast example;
#' # a real analysis would use the full cohort.
#' sample_ids <- setdiff(colnames(mrna_data$groupA), "gene_name")[1:10]
#' mrna_sample_data <- mrna_data$groupA[, c("gene_name", sample_ids)]
#' protein_sample_data <- protein_data$groupA[, c("ref_seq", "gene_name", sample_ids)]
#'
#' # A few genes have missing measurements for some samples in this subset;
#' # they are handled like in the rest of DrDimont, via handling_missing_data.
#' example_settings <- DrDimont::drdimont_settings(
#'   handling_missing_data = "pairwise.complete.obs"
#' )
#'
#' \donttest{
#' result <- create_lioness_multiomics_samples(
#'   sample_layers = list(
#'     mrna = list(raw = mrna_sample_data, meta_cols = "gene_name"),
#'     protein = list(raw = protein_sample_data, meta_cols = c("ref_seq", "gene_name"))
#'   ),
#'   # Lowered from the default (0.8): DrDimont's demo panel has only
#'   # ~50 genes, too few for a strict scale-free fit on every single-sample
#'   # network. Use a stricter cutoff for real, genome-scale data.
#'   r_squared_cutoff_samples = 0.1,
#'   settings = example_settings,
#'   connect_on = "gene_name"
#' )
#' }
#'
#' \dontrun{
#' # Comparing two samples' drug response additionally requires drug-target
#' # interactions and a configured Python environment for the interaction-score
#' # step (see DrDimont::install_python_dependencies()), like the rest of
#' # DrDimont's pipeline.
#' data(drug_gene_interactions, package = "DrDimont")
#' drug_target_input <- DrDimont::make_drug_target(
#'   target_molecules = "protein",
#'   interaction_table = drug_gene_interactions,
#'   match_on = "gene_name"
#' )
#'
#' result <- create_lioness_multiomics_samples(
#'   sample_layers = list(
#'     mrna = list(raw = mrna_sample_data, meta_cols = "gene_name"),
#'     protein = list(raw = protein_sample_data, meta_cols = c("ref_seq", "gene_name"))
#'   ),
#'   drug_target_interactions = drug_target_input,
#'   compare = list(
#'     type = "sample_vs_sample",
#'     samples = list(c(sample_ids[1], sample_ids[2]))
#'   ),
#'   r_squared_cutoff_samples = 0.1,
#'   settings = example_settings,
#'   connect_on = "gene_name"
#' )
#' }
#'
#' @export
create_lioness_multiomics_samples <- function(
    sample_layers,
    group_layers = NULL,
    drug_target_interactions = NULL,
    compare = NULL,
    r_squared_cutoff_samples = 0.8,
    r_squared_cutoff_groups = 0.8,
    settings = DrDimont::drdimont_settings(),
    connect_on = NULL,
    clip_LIONESS_weights = TRUE,
    min_layer = 2,
    min_overlap = 0.7,
    save_intermediate = FALSE,
    intermediate_dir = NULL,
    save_results = FALSE,
    results_dir = NULL,
    group_reduction_method = "pickHardThreshold"
) {
  # The skeletal structure of these comments was based on the comments structure of DrDimont, explained by ChatGPT (LLM1)

  # Checks sample_layers, group_layers and connect_on the same way DrDimont checks its input (check_lioness_input in checks.R).
  # All problems are collected and reported together, then the workflow stops.
  DrDimont::return_errors(check_lioness_input(sample_layers, group_layers, connect_on))
  
  # Checks if there was a intermediate_dir value passed in the function call, if the user wants to save intermediates, stops and warns if not the case.
  # If directory has been passed it is created and it is checked if it does now exist, if not stop and warning.
  if (isTRUE(save_intermediate)) {
    if (is.null(intermediate_dir)) {
      stop("intermediate_dir must be provided when save_intermediate = TRUE.", call. = FALSE)
    }
    dir_created <- dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
    if (!isTRUE(dir_created) && (!dir.exists(intermediate_dir))) {
      stop(paste0("intermediate_dir does not exist and could not be created: ", intermediate_dir), call. = FALSE)
    }
  }
  
  # Checks if there was a results_dir value passed in the function call, if the user wants to save results, stops and warns if not the case.
  # If directory has been passed it is created and it is checked if it does now exist, if not stop and warning.
  if (isTRUE(save_results)) {
    if (is.null(results_dir)) {
      stop("results_dir must be provided when save_results = TRUE.", call. = FALSE)
    }
    dir_created <- dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
    if (!isTRUE(dir_created) && (!dir.exists(results_dir))) {
      stop(paste0("results_dir does not exist and could not be created: ", results_dir), call. = FALSE)
    }
  }
  
  # First checks if the user wants to compare, if yes, and if at least one group vs sample comparison is present:
  # Checks if group and sample data have the same genes in the same order because later node_ids are generated for the graphs of the networks and they need to be in the same order as well.
  # Identifiers containing "_" are dropped because LIONESS interaction names are later split by this character. Duplicates are also dropped.
  if (!is.null(compare) && any(compare$type == "sample_vs_group")) {
    aligned <- align_sample_group_raw_layers(
      sample_layers = sample_layers,
      group_layers = group_layers,
      connect_on = connect_on,
      min_overlap = min_overlap
    )
    sample_layers <- aligned$sample_layers
    group_layers  <- aligned$group_layers
  }
  
  # The raw data files get transformed into a form LIONESS recognizes, all column with non numerical data (all the meta_cols) get removed.
  # At the end only columns with numerical data and the unique identifier column (name_col) remain (with the entries in the name_col as row names).
  # Identifiers containing "_" are dropped because LIONESS interaction names are later split by this character. Duplicates are also dropped.
  expr_data <- lapply(sample_layers, function(layer) {
    prepare_lioness_expr(
      raw_data = layer$raw,
      identifier_col = connect_on,
      meta_cols = layer$meta_cols
    )
  })

  # Find sample columns for each layer.
  sample_names_by_layer <- lapply(expr_data, colnames)
  
  # Count in how many layers each sample occurs.
  sample_counts <- table(unlist(sample_names_by_layer, use.names = FALSE))
  
  # Save names of samples that occur in at least two layers in variable.
  common_samples <- names(sample_counts)[sample_counts >= min_layer]
  
  # stop and warn user if no samples occur in at least min_layer.
  if (length(common_samples) == 0) {
    stop(paste0("No samples were found in at least ", min_layer, " layers."), call. = FALSE)
  }
  
  #LIONESS is run to create a single sample molecular network for all samples, for this the default aggregate network upon which LIONESS builds its single sample molecular networks is PEARSON.
  # Missing values are handled per layer with the handling_missing_data setting, the same way DrDimont handles them for the group networks.
  lioness_samples <- lapply(names(expr_data), function(layer_name) {
    run_lioness(
      expr = expr_data[[layer_name]],
      handling_missing_data = get_layer_setting(layer_name, "", settings, "handling_missing_data")
    )
  })
  names(lioness_samples) <- names(expr_data)
  
  # Converts the sample list from LIONESS (with interaction (GeneA to GeneB) in one column and weight in another) to a adjacency matrix with the Genes in row1 and column1.
  # Also removes edges bigger than 1 and smaller than -1 or scales them to a range between -1 and 1, because DrDimont can only work in this range.
  adj_matrices <- lapply(lioness_samples, function(sample_list) {
    sample_list <- sample_list[names(sample_list) %in% common_samples]
    lapply(sample_list, function(sample) {
      edge_to_adjacency(
        sample = sample,
        clip_LIONESS_weights = clip_LIONESS_weights
      )
    })
  })
  
  # Reduces the networks via the DrDimont method pickHardThreshold.
  reduced_samples <- lapply(adj_matrices, function(sample_adj_list) {
    lapply(sample_adj_list, function(adj_mat) {
      result <- network_reduction_by_pickHardThreshold(
        adjacency_matrix = adj_mat,
        r_squared_cutoff = r_squared_cutoff_samples,
        cut_vector = seq(0.01, 0.99, by = 0.01)
      )
      # DrDimont's reduction marks removed edges with -999, they are set to 0 so they are left out of the graph.
      result[result == -999] <- 0
      return(result)
    })
  })
  
  # Builds multi-omics networks with samples that exist in at least two layers (common_samples), this functionality was optimized by AI (LLM3).
  multiomics_samples <- lapply(common_samples, function(sample_id) {
    
    # Keep only layers where this sample exists and save them in sample_adjs.
    sample_adjs <- reduced_samples[vapply(reduced_samples, function(layer_samples) {sample_id %in% names(layer_samples)}, logical(1))]
    
    # Now every layer gets reduced to only hold the sample_id and not every other sample as well.
    sample_adjs <- lapply(sample_adjs, function(layer_samples) {
      layer_samples[[sample_id]]
    })
    
    # multi-omics network gets build for specific sample.
    build_sample_multilayer(
      sample_adjs = sample_adjs,
      settings = settings,
      identifier_col = connect_on
    )
  })
  
  # The samples get their specific names from the common_samples variable.
  names(multiomics_samples) <- common_samples
  
  # If intermediate data should be saved, the new multi-omics networks are saved as .rds object.
  if(isTRUE(save_intermediate)){
    result_saver(
      object = multiomics_samples,
      filename = "multiomics_samples.rds",
      dir = intermediate_dir
    )
  }
  
  # Define comparison_result, interaction score cache and group multiomics, because if the check below is negative it would be undefined.
  comparison_result <- NULL
  group_multiomics <- NULL
  interaction_score_cache <- NULL
  
  # Checks if there were variables passed in the function call to compare, skips the compare step if none were passed.
  if (!is.null(compare)) {
    
    # Checks if drug target interactions were passed.
    if (is.null(drug_target_interactions)) {
      stop("drug_target_interactions must be provided when compare is used.", call. = FALSE)
    }
    
    # Check what is to be compared to create group network if needed.
    if (any(compare$type == "sample_vs_group")) {
      
      # Check if group layers were passed.
      if (is.null(group_layers)) {
        stop("group_layers must be provided for sample_vs_group comparison.", call. = FALSE)
      }
      
      # Set group settings
      group_settings <- settings
      
      # Set reduction method for group, depending on what method was chosen, default is pickHardThreshold and throw an error and stop if wrong value was input.
      if (group_reduction_method == "pickHardThreshold") {
        group_settings$r_squared_cutoff <- r_squared_cutoff_groups
        group_settings$cut_vector <- seq(0.01, 0.99, by = 0.01)
      } else if (group_reduction_method == "p_value") {
        group_settings$reduction_method <- "p_value"
        group_settings$p_value_adjustment_method <- "BH"
        group_settings$reduction_alpha <- 0.05
      } else{
        stop("group_reduction_method can only be set to 'pickHardThreshold' or 'p_value'.", call. = FALSE)
      }
      
      # Build a multi-omics network for each group and save them in the variable group_multiomics.
      for (group_name in names(group_layers)) {
        group_multiomics[[group_name]] <- build_group_multilayer_drdimont(
          group_layers = group_layers[[group_name]],
          settings = group_settings,
          connect_on = connect_on
        )
      }
      
      #If the user want to save intermediate results, the group multi-omics networks are saved.
      if(isTRUE(save_intermediate)){
        result_saver(
          object = group_multiomics,
          filename = "multiomics_group.rds",
          dir = intermediate_dir
        )
      }
    }
    
    # Create variable to store interaction scores.
    interaction_score_cache <- list()
    
    # Find out which networks are to be compared and save them in variable unique_net_ids, use.names has been set to FALSE to avoid labels.
    # The networks that are to be compared are only the unique ones, so the interaction scores that are computed are only computed once per network.
    unique_net_ids <- unique(unlist(compare$samples, use.names = FALSE))
    
    # Tell user interaction scores are computed right now.
    message("Pre-computing interaction scores for ", length(unique_net_ids), " unique network(s): ", paste0(unique_net_ids, collapse = ", "))
    
    # Iterate through every network, create interaction scores for them and save them in interaction_score_cache.
    for (net_id in unique_net_ids) {
      message("  Computing scores for: ", net_id)
      
      # Check if network is a group network or sample network, and if it is a sample network if it occurs in common_samples, stop if not found there, because it did not occur in enough layers then or user has made a typo.
      # If networks exist, save the multi-omics networks of them in variable network.
      if (net_id %in% names(group_multiomics)) {
        network <- group_multiomics[[net_id]]
      } else {
        if (!net_id %in% common_samples) {
          stop(paste0("Sample or group '", net_id, "' was not found in at least ", min_layer, " layers."), call. = FALSE)
        }
        network <- multiomics_samples[[net_id]]
      }
      
      # Calculate interactions score graphs for every unique multi-omics network to be compared and save them in variable interaction_score_cache.
      interaction_score_cache[[net_id]] <- compute_single_network_scores(
        network = network,
        drug_target_interactions = drug_target_interactions,
        settings = settings
      )
    }
    
    # Tell user interaction scores have been computed and comparisons have started.
    message("Pre-computation complete. Running comparisons...")
    
    # Compare everything to be compared
    for (i in seq_along(compare$type)) {
      
      # Declare type of comparison (comparison_type) and samples or groups to be compared (comparison_samples).
      comparison_type    <- compare$type[i]
      comparison_samples <- compare$samples[[i]]
      
      # If the User wants to compare two samples this code is executed.
      # Otherwise the else if statement holds true or an error is thrown and the code is stopped, for example if a typo happens.
      if (comparison_type == "sample_vs_sample") {
        
        # Check if two samples were provided for comparison, if not stop and throw an error.
        if (length(comparison_samples) != 2) {
          stop("For sample_vs_sample, exactly two samples must be provided.", call. = FALSE)
        }
        
        # Declare the two samples to be compared and put them in the variables sample_a and sample_b.
        sample_a <- comparison_samples[1]
        sample_b <- comparison_samples[2]
        
        # Create a result_name variable, with which the specific comparison can be found later on.
        result_name <- paste0("sample_vs_sample__ ", sample_a, " __vs__ ", sample_b)
        
        # Run the comparison and save the result with the result_name as name in the comparison_result variable.
        comparison_result[[result_name]] <- compare_drug_interactions(
          drug_target_interactions = drug_target_interactions,
          settings = settings,
          precomputed_a = interaction_score_cache[[sample_a]],
          precomputed_b = interaction_score_cache[[sample_b]]
        )
        
      } else if (comparison_type == "sample_vs_group") {
        
        # Check if two inputs were provided for comparison, if not stop and throw an error.
        if (length(comparison_samples) != 2) {
          stop("For sample_vs_group, exactly one sample and one group must be provided.", call. = FALSE)
        }

        # Declare the sample and group to be compared and put them in the variables sample and group.
        sample <- comparison_samples[1]
        group <- comparison_samples[2]
        
        # Validate that the connect_on column is in the same order between the sample and group annotations (required because node_ids are position-based).
        sample_ann <- multiomics_samples[[sample]]$annotations$both
        group_ann  <- group_multiomics[[group]]$annotations$both
        for (ln in unique(sample_ann$layer)) {
          validate_connect_on_alignment(sample_ann, group_ann, connect_on, ln)
        }
        
        # Create a result_name variable, with which the specific comparison can be found later on.
        result_name <- paste0("sample_vs_group__ ", sample, " __vs__ ", group)
        
        # Run the comparison and save the result with the result_name as name in the comparison_result variable.
        comparison_result[[result_name]] <- compare_drug_interactions(
          drug_target_interactions = drug_target_interactions,
          settings = settings,
          precomputed_a = interaction_score_cache[[sample]],
          precomputed_b = interaction_score_cache[[group]]
        )
      } else{
        stop(paste0("comparison_type must either be sample_vs_group or sample_vs_sample, but '", comparison_type, "' was passed"), call. = FALSE)
      }
    }
  }
  
  # Save the results in the variable final_result.
  final_result <- list(
    multiomics_per_sample = multiomics_samples,
    comparison_result = comparison_result,
    group_multiomics = group_multiomics,
    interaction_score_cache = interaction_score_cache
  )
  
  # Save the final_result as a .rds object if User wants to save the final output.
  if (isTRUE(save_results)) {
    result_saver(
      object = final_result,
      filename = "lioness_multiomics_result.rds",
      dir = results_dir
    )
  }
  return(final_result)
}