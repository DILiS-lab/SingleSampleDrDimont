# The two functions in this file are copied unchanged from DrDimont 0.1.7 (CRAN,
# published 2026-07-21), get_layer_setting() from R/settings.R and network_reduction_by_pickHardThreshold()
# from R/network_reduction.R. DrDimont does not export them, and CRAN does not allow
# calling unexported functions of another package with `:::`.
# DrDimont is published under the MIT license, Copyright (c) 2026 DrDimont authors
# (see LICENSE.md).

get_layer_setting <- function(layer, group, settings, setting_name) {
    #' @title [INTERNAL] Get layer (and group) settings
    #'
    #' @description Returns specified setting for a specific network layer (and group).
    #'
    #' @param layer [list] A network layer created by \code{\link[DrDimont]{make_layer}}
    #' @param group [string] A network group
    #' @param settings [list] Named list of settings created by \code{\link[DrDimont]{drdimont_settings}}
    #' @param setting_name [string] String indicating the setting to return.
    #' @return Setting value(s) for this layer (and group)
    #' 
    #' @keywords internal
    #' @noRd
    
    if (!is.list(settings[[setting_name]])) {return(settings[[setting_name]])}
    else if (!is.null(settings[[setting_name]][[group]][[layer]])) {return(settings[[setting_name]][[group]][[layer]])}
    else if (!is.null(settings[[setting_name]][[layer]])) {return(settings[[setting_name]][[layer]])}
    else if (!is.null(settings[[setting_name]][[group]][['default']])) {return(settings[[setting_name]][[group]][['default']])}
    else if (!is.null(settings[[setting_name]][['default']])){return(settings[[setting_name]][['default']])}
    else {stop(stringr::str_interp("Neither was a setting ${setting_name} given for layer ${layer} nor any default."))}
}


network_reduction_by_pickHardThreshold <- function(adjacency_matrix,
                                                   r_squared_cutoff=0.85,
                                                   cut_vector=seq(0.2, 0.8, by = 0.01),
                                                   mean_number_edges=NULL,
                                                   edge_density=NULL) {
    #' @title [INTERNAL] Reduces network based on WGCNA::pickHardThreshold function
    #'
    #' @description  [INTERNAL] This function uses \code{\link[WGCNA]{pickHardThreshold.fromSimilarity}} to analyze
    #' scale free topology for multiple hard thresholds. A cutoff is estimated, if no cutoff is
    #' found the function terminates with an error message. All values below the cutoff will be set to NA and the
    #' reduced adjacency is returned.
    #'
    #' @param adjacency_matrix [matrix] Adjacency matrix of correlations computed using \code{\link[WGCNA]{cor}} in
    #' \code{\link[DrDimont]{compute_correlation_matrices}}
    #' @param r_squared_cutoff [float] A number indicating the desired minimum scale free topology fitting index R^2 for reduction
    #' using \code{\link[WGCNA]{pickHardThreshold}}. (default: 0.85)
    #' @param cut_vector [sequence of float] A vector of hard threshold cuts for which the scale free topology fit indices are to
    #' be calculated during reduction with \code{\link[WGCNA]{pickHardThreshold}}. (default: seq(0.2, 0.8, by = 0.01))
    #' @param mean_number_edges [int] Find a suitable edge weight cutoff employing \code{\link[WGCNA]{pickHardThreshold}} to reduce
    #' the network to at most the specified mean number of edges. Attention: This parameter overwrites the 'r_squared_cutoff' and
    #' 'edge_density' parameters if not set to NULL. (default: NULL)
    #' @param edge_density [float] Find a suitable edge weight cutoff employing \code{\link[WGCNA]{pickHardThreshold}} to reduce the
    #' network to at most the specified edge density. Attention: This parameter overwrites the 'r_squared_cutoff' parameter if not set
    #' to NULL. (default: NULL)
    #' 
    #' @source The original implementation of pickHardThreshold is used from
    #' \code{\link[WGCNA]{pickHardThreshold.fromSimilarity}}
    #'
    #' @return A reduced adjacency matrix of correlations with NA's inserted at positions below
    #' estimated cutoff.
    #' 
    #' @keywords internal
    #' @noRd
    
    message(format(Sys.time(), "[%y-%m-%d %X] "), 'Reducing network by WGCNA::pickHardThreshold...')

    ### if mean number of edges given calculate cut threshold based on WGCNA::pickHardThreshold.fromSimilarity() mean.k. values
    ### else use the r_squared_cutoff
    if (!is.null(mean_number_edges)){

        message(format(Sys.time(), "[%y-%m-%d %X] "), 'Mean number of edges: ', as.character(mean_number_edges))

        n_entities <- nrow(adjacency_matrix)

        ### "invisible(capture.output())" suppresses printing of a table of cutoff and R^2 values
        ### WGCNA pickHardThreshold computation
        invisible(utils::capture.output(
            wgcna_tresholds <- WGCNA::pickHardThreshold.fromSimilarity(abs(adjacency_matrix), cutVector=cut_vector)
        ))

        ### get WGCNA cut with mean number of edges <= mean_number_edges
        cut_estimate = NULL
        for (row in 1:nrow(wgcna_tresholds$fitIndices)) {
            cut_val <- wgcna_tresholds$fitIndices[row, "Cut"]
            r2_val  <- wgcna_tresholds$fitIndices[row, "SFT.R.sq"]
            slope_val  <- wgcna_tresholds$fitIndices[row, "slope."]
            mean_k_val <- wgcna_tresholds$fitIndices[row, "mean.k."]

            if ((mean_k_val <= mean_number_edges*2)){ #(slope_val<0) & 
                cut_estimate=cut_val
                message(format(Sys.time(), "[%y-%m-%d %X] "), 'R2 cutoff: ', round(r2_val, 2))
                message(format(Sys.time(), "[%y-%m-%d %X] "), 'Cut Threshold: ', as.character(cut_estimate))
                break
            }
        }

        ### terminate execution if pickHardThreshold cannot find cutoff
        if (is.null(cut_estimate)) {
            message(format(Sys.time(), "[%y-%m-%d %X] "),
                    "ERROR: WGCNA::pickHardThreshold: failed to find a threshold with given mean number of edges per node and cut vector.
                    Please, try a different mean number of edges or cut vector. Lowest mean number of edges computed was ",
                    round(min(wgcna_tresholds$fitIndices[wgcna_tresholds$fitIndices$slope.<0, "mean.k."], na.rm=T)/2, 3),
                    ". Highest mean number of edges computed was ",
                    round(max(wgcna_tresholds$fitIndices[wgcna_tresholds$fitIndices$slope.<0, "mean.k."], na.rm=T)/2, 3), ".")
            stop("WGCNA::pickHardThreshold: failed to find a threshold with given mean number of edges per node and cut vector.
                 Please, try a different mean number of edges or cut vector. Lowest mean number of edges computed was ",
                 round(min(wgcna_tresholds$fitIndices[wgcna_tresholds$fitIndices$slope.<0, "mean.k."], na.rm=T)/2, 3),
                 ". Highest mean number of edges computed was ",
                 round(max(wgcna_tresholds$fitIndices[wgcna_tresholds$fitIndices$slope.<0, "mean.k."], na.rm=T)/2, 3), ".")
        }

    ### if mean number of edges not given and edge density given calculate cut threshold based on WGCNA::pickHardThreshold.fromSimilarity() mean.k. values
    } else if (!is.null(edge_density)){

        message(format(Sys.time(), "[%y-%m-%d %X] "), 'Network edge density: ', as.character(edge_density))

        n_entities <- nrow(adjacency_matrix)

        ### "invisible(capture.output())" suppresses printing of a table of cutoff and R^2 values
        ### WGCNA pickHardThreshold computation
        invisible(utils::capture.output(
            wgcna_tresholds <- WGCNA::pickHardThreshold.fromSimilarity(abs(adjacency_matrix), cutVector=cut_vector)
        ))

        ### get WGCNA cut with edge density <= edge_density
        cut_estimate <- NULL
        for (row in 1:nrow(wgcna_tresholds$fitIndices)) {
            cut_val <- wgcna_tresholds$fitIndices[row, "Cut"]
            r2_val  <- wgcna_tresholds$fitIndices[row, "SFT.R.sq"]
            slope_val  <- wgcna_tresholds$fitIndices[row, "slope."]
            mean_k_val <- wgcna_tresholds$fitIndices[row, "mean.k."]

            if ((slope_val<0) & (mean_k_val <= (edge_density*n_entities-edge_density)*2)){
                cut_estimate=cut_val
                message(format(Sys.time(), "[%y-%m-%d %X] "), 'R2 cutoff: ', round(r2_val, 2))
                message(format(Sys.time(), "[%y-%m-%d %X] "), 'Cut Threshold: ', as.character(cut_estimate))
                break
            }
        }

        ### terminate execution if pickHardThreshold cannot find cutoff
        if (is.null(cut_estimate)) {
            message(format(Sys.time(), "[%y-%m-%d %X] "),
                    "ERROR: WGCNA::pickHardThreshold: failed to find a threshold with given edge density and cut vector.
                    Please, try a different edge density or cut vector. Lowest edge density computed was ",
                    round(min(wgcna_tresholds$fitIndices[wgcna_tresholds$fitIndices$slope.<0, "mean.k."], na.rm=T)/(2*(n_entities-1)), 5),
                    ". Highest edge density computed was ",
                    round(max(wgcna_tresholds$fitIndices[wgcna_tresholds$fitIndices$slope.<0, "mean.k."], na.rm=T)/(2*(n_entities-1)), 3), ".")
            stop("WGCNA::pickHardThreshold: failed to find a threshold with given edge density and cut vector.
                 Please, try a different edge density or cut vector. Lowest edge density computed was ",
                 round(min(wgcna_tresholds$fitIndices[wgcna_tresholds$fitIndices$slope.<0, "mean.k."], na.rm=T)/(2*(n_entities-1)), 5),
                 ". Highest edge density computed was ",
                 round(max(wgcna_tresholds$fitIndices[wgcna_tresholds$fitIndices$slope.<0, "mean.k."], na.rm=T)/(2*(n_entities-1)), 3), ".")
        }
    ### else use the r_squared_cutoff
    } else {

        message(format(Sys.time(), "[%y-%m-%d %X] "), 'R2 cutoff: ', as.character(r_squared_cutoff))
        ### "invisible(capture.output())" suppresses printing of a table of cutoff and R^2 values
        ### WGCNA pickHardThreshold computation
        invisible(utils::capture.output(
            wgcna_tresholds <- WGCNA::pickHardThreshold.fromSimilarity(abs(adjacency_matrix), r_squared_cutoff, cut_vector)
        ))

        cut_estimate = wgcna_tresholds$cutEstimate
        ### terminate execution if pickHardThreshold cannot find cutoff
        if (is.na(cut_estimate)) {
            message(format(Sys.time(), "[%y-%m-%d %X] "),
                    "ERROR: WGCNA::pickHardThreshold failed to find a suitable cutoff with the given cut_vector at the given R^2 cutoff.
                    Please, try a different cut_vetor or a different cutoff. Highest R^2 cutoff computed was ",
                    round(max(wgcna_tresholds$fitIndices[wgcna_tresholds$fitIndices$slope.<0, "SFT.R.sq"], na.rm=T), 3), ".")
            stop("WGCNA::pickHardThreshold failed to find a suitable cutoff with the given cut_vector at the given R^2 cutoff.
                    Please, try a different cut_vetor or a different cutoff. Highest R^2 cutoff computed was ",
                 round(max(wgcna_tresholds$fitIndices[wgcna_tresholds$fitIndices$slope.<0, "SFT.R.sq"], na.rm=T), 3), ".")
        ##### check if cut estimate was found and if slope at cut is negative
        ##### if true continue, else find a cut estimate with negative slope
        } else if ((is.na(cut_estimate)==FALSE) & (wgcna_tresholds$fitIndices[wgcna_tresholds$fitIndices$Cut==cut_estimate,]$slope < 0)){
            message(format(Sys.time(), "[%y-%m-%d %X] "), 'Cut Threshold: ', as.character(cut_estimate))
        } else {
            cut_estimate=NA
            for (row in 1:nrow(wgcna_tresholds$fitIndices)) {
                cut_val <- wgcna_tresholds$fitIndices[row, "Cut"]
                r2_val  <- wgcna_tresholds$fitIndices[row, "SFT.R.sq"]
                slope_val  <- wgcna_tresholds$fitIndices[row, "slope."]

                if ((r2_val>r_squared_cutoff) & (slope_val<0)){
                    cut_estimate=cut_val
                    message(format(Sys.time(), "[%y-%m-%d %X] "), 'Cut Threshold: ', as.character(cut_estimate))
                    break
                }
            }
        }

        ### terminate execution if pickHardThreshold cannot find cutoff
        if(is.na(cut_estimate)){ #Please try setting a lower R-squared or using a different network reduction method. Highest R-squared computed was ",round(highestR2, 3)
            message(format(Sys.time(), "[%y-%m-%d %X] "),
                    "ERROR: WGCNA::pickHardThreshold failed to find a suitable cutoff with the given cut_vector at the given R^2 cutoff.
                    Please, try a different cut_vetor or a different cutoff. Highest R^2 cutoff computed was ",
                    round(max(wgcna_tresholds$fitIndices[wgcna_tresholds$fitIndices$slope.<0, "SFT.R.sq"], na.rm=T), 3), ".")
            stop("WGCNA::pickHardThreshold failed to find a suitable cutoff with the given cut_vector at the given R^2 cutoff.
                    Please, try a different cut_vetor or a different cutoff. Highest R^2 cutoff computed was ",
                 round(max(wgcna_tresholds$fitIndices[wgcna_tresholds$fitIndices$slope.<0, "SFT.R.sq"], na.rm=T), 3), ".")
        }
    }

    cutEstimate <- cut_estimate

    ### set all entries below the threshold to NA
    adjacency_matrix[abs(adjacency_matrix) < cutEstimate] <- -999

    return(adjacency_matrix)
}
