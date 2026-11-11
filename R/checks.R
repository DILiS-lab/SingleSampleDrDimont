check_lioness_input <- function(sample_layers, group_layers, connect_on) {
    #' @title [INTERNAL] Check LIONESS workflow input
    #'
    #' @description [INTERNAL] Checks if the layers and the identifier column given to
    #' \code{\link{create_lioness_multiomics_samples}} are valid and have the right format
    #'
    #' @param sample_layers [list] Named list of sample layers, each containing `raw` and `meta_cols`
    #' @param group_layers [list|NULL] Named list of groups, each a named list of layers like `sample_layers`
    #' @param connect_on [string] Name of the identifier column
    #' @return Character string vector containing error messages.
    #'
    #' @keywords internal
    #' @noRd

    errors <- c()

    ### check if argument 'connect_on' is a character string of length 1
    if (!is.character(connect_on) || length(connect_on) != 1) {
        errors <- c(errors, "The argument `connect_on` has to be declared as one character string naming the identifier column, e.g. \"gene_name\".")
        connect_on <- NULL
    }

    ### collect all layers to check, labelled by where they come from
    layers <- list()
    if (!is.list(sample_layers) || length(sample_layers) == 0 || is.null(names(sample_layers)) || any(names(sample_layers) == "")) {
        errors <- c(errors, "The argument `sample_layers` has to be a named list of layers, e.g. list(mrna = list(raw = ..., meta_cols = ...)).")
    } else {
        for (layer_name in names(sample_layers)) {
            layers[[stringr::str_interp("Sample layer \"${layer_name}\"")]] <- sample_layers[[layer_name]]
        }
    }
    if (!is.null(group_layers)) {
        if (!is.list(group_layers) || length(group_layers) == 0 || is.null(names(group_layers)) || any(names(group_layers) == "")) {
            errors <- c(errors, "The argument `group_layers` has to be a named list of groups, each a named list of layers, e.g. list(group1 = list(mrna = list(raw = ..., meta_cols = ...))).")
        } else {
            for (group_name in names(group_layers)) {
                group <- group_layers[[group_name]]
                if (!is.list(group) || length(group) == 0 || is.null(names(group)) || any(names(group) == "")) {
                    errors <- c(errors, stringr::str_interp("Group \"${group_name}\": The group has to be a named list of layers, e.g. list(mrna = list(raw = ..., meta_cols = ...))."))
                    next
                }
                for (layer_name in names(group)) {
                    layers[[stringr::str_interp("Group \"${group_name}\", layer \"${layer_name}\"")]] <- group[[layer_name]]
                }
            }
        }
    }

    for (label in names(layers)) {
        layer <- layers[[label]]

        ### check if layer contains 'raw' and 'meta_cols'
        if (!is.list(layer) || !all(c("raw", "meta_cols") %in% names(layer))) {
            errors <- c(errors, stringr::str_interp("${label}: The layer has to be a list containing `raw` and `meta_cols`."))
            next
        }

        ### check if 'raw' is a dataframe
        if (!is.data.frame(layer$raw)) {
            errors <- c(errors, stringr::str_interp("${label}: `raw` has to be a dataframe."))
            next
        }
        column_names <- colnames(layer$raw)

        ### check if 'meta_cols' is a character vector naming columns of 'raw'
        if (!is.character(layer$meta_cols)) {
            errors <- c(errors, stringr::str_interp("${label}: `meta_cols` has to be a character vector naming the non numerical columns of `raw`."))
            next
        }
        missing_meta_cols <- setdiff(layer$meta_cols, column_names)
        if (length(missing_meta_cols) > 0) {
            errors <- c(errors, stringr::str_interp("${label}: The columns in `meta_cols` were not found in `raw`: ${toString(missing_meta_cols)}. Column names are: ${toString(column_names)}."))
        }

        ### check if the identifier column given in 'connect_on' is a metadata column of 'raw'
        if (!is.null(connect_on)) {
            if (!connect_on %in% column_names) {
                errors <- c(errors, stringr::str_interp("${label}: The identifier column specified in `connect_on` was not found in `raw`. Argument `connect_on` is set to: ${connect_on}. Column names are: ${toString(column_names)}."))
            } else if (!connect_on %in% layer$meta_cols) {
                errors <- c(errors, stringr::str_interp("${label}: The identifier column specified in `connect_on` has to be listed in `meta_cols`. Argument `connect_on` is set to: ${connect_on}."))
            }
        }

        ### check if the data without the metadata columns is numeric (same check as in check_layer)
        if (!is.numeric(as.matrix(layer$raw[, setdiff(column_names, layer$meta_cols), drop = FALSE]))) {
            errors <- c(errors, stringr::str_interp("${label}: The supplied data is not numeric. Please check your data. Numerical data is required."))
        }
    }

    return(errors)
}
