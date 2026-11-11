# SingleSampleDrDimont


## Description

Molecular networks built from omics data, such as mRNA or protein expression, describe how
molecules interact and can be used to study diseases. They are usually built for groups of
patients or samples to describe a condition, which can hide the properties of an individual
patient or sample behind the group average. The package enables individual-sample comparisons.
It builds single-sample networks for one or more omics layers with the LIONESS method (Linear
Interpolation to Obtain Network Estimates for Single Samples) of Kuijjer et al. (2019)
<https://doi.org/10.1016/j.isci.2019.03.021>, as implemented in 'lionessR', and integrates them
into multi-omics networks with 'DrDimont' by Hiort et al. (2022)
<https://doi.org/10.1093/bioinformatics/btac477>. The single-sample multi-omics networks can be
compared with each other or with a group network to predict differences in drug response for
individual samples or patients.

## Installation

### Installing the R package

Four dependencies come from Bioconductor: `lionessR` and
`SummarizedExperiment`, plus `impute` and `preprocessCore`, which `WGCNA` needs.
As `install.packages()` only searches CRAN, install with `BiocManager`.

- From CRAN:
  ```r
  install.packages("BiocManager")
  BiocManager::install("SingleSampleDrDimont")
  ```
- From source: either clone the repository and use `devtools::install()` within R, or use
  `remotes::install_github("DILiS-lab/SingleSampleDrDimont")` without cloning. Both install the Bioconductor
  dependencies automatically.

### Installing the Python dependencies
To compare samples, a Python (tested with Python 3.10) installation is required. 
Afterwards, use `DrDimont::install_python_dependencies()` to install the
necessary dependencies. The `package_manager` argument selects `"pip"` (
default) or `"conda"`. If you use conda, also set `conda = TRUE` in
`DrDimont::drdimont_settings()`.

If you prefer to install the dependencies manually, the requirements files ship with DrDimont:
`system.file("requirements_pip.txt", package = "DrDimont")` and
`system.file("requirements_conda.txt", package = "DrDimont")`. Install them into a virtualenv or
conda environment named `r-DrDimont` and use python-igraph>= 1.0.0, which the automatic installer also requires.

## Example

Build single-sample multi-omics networks for 10 samples of the example data that ships
with DrDimont:

```r
library(SingleSampleDrDimont)

data(mrna_data, package = "DrDimont")
data(protein_data, package = "DrDimont")

sample_ids <- setdiff(colnames(mrna_data$groupA), "gene_name")[1:10]
mrna_sample_data <- mrna_data$groupA[, c("gene_name", sample_ids)]
protein_sample_data <- protein_data$groupA[, c("ref_seq", "gene_name", sample_ids)]

result <- create_lioness_multiomics_samples(
  sample_layers = list(
    mrna = list(raw = mrna_sample_data, meta_cols = "gene_name"),
    protein = list(raw = protein_sample_data, meta_cols = c("ref_seq", "gene_name"))
  ),
  r_squared_cutoff_samples = 0.1,
  settings = DrDimont::drdimont_settings(handling_missing_data = "pairwise.complete.obs"),
  connect_on = "gene_name"
)
```

See `?create_lioness_multiomics_samples` for the comparison options.

## References

This package builds on:
- DrDimont: Hiort et al. (2022), Bioinformatics 38(Supplement_2), ii113–ii119,
  https://doi.org/10.1093/bioinformatics/btac477
- LIONESS method: Kuijjer et al. (2019), iScience 14, 226–240,
  https://doi.org/10.1016/j.isci.2019.03.021
- lionessR: Kuijjer et al. (2019), BMC Cancer,
  https://doi.org/10.1186/s12885-019-6235-7

## License

MIT, see `LICENSE.md`.
