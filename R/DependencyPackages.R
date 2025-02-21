# Install 'devtools' if it's not already installed
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")

library(devtools)

# Install 'BiocManager' if it's not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Function to install missing CRAN packages
install_if_missing <- function(package_list) {
  missing_packages <- package_list[!sapply(package_list, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    install.packages(missing_packages)
  }
}

# Function to install missing Bioconductor packages
install_bioc_if_missing <- function(package_list) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  missing_packages <- package_list[!sapply(package_list, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    BiocManager::install(missing_packages)
  }
}

# Function to install missing GitHub packages
install_github_if_missing <- function(repo, ref = NULL, ...) {
  package_name <- unlist(strsplit(repo, "/"))[2]  # Extract package name from repo string
  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
  
  if (!requireNamespace(package_name, quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE))
      install.packages("devtools")
    
    if (!is.null(ref)) {
      devtools::install_github(repo, ref = ref, ...)
    } else {
      devtools::install_github(repo, ...)
    }
  }
}

# List of required CRAN packages
required_cran_packages <- c("Seurat", "dplyr", "ggplot2", "shadowtext", "scales", "cowplot",
                            "data.table", "Matrix", "matrixStats", "fs", "gridExtra",
                            "magrittr", "tibble", "lsa")

# List of required Bioconductor packages
required_bioc_packages <- c("SingleCellExperiment", "HDF5Array", "SpatialExperiment",
                            "bluster", "BiocParallel", "BioQC", "NMI", "coop", "scater")

# Install CRAN packages if they are missing
install_if_missing(required_cran_packages)

# Install Bioconductor packages if they are missing
install_bioc_if_missing(required_bioc_packages)

# Install `scater` with all dependencies
if (!requireNamespace("scater", quietly = TRUE)) {
  BiocManager::install("scater", dependencies = TRUE)
}

#  Only check SpatialFeatureExperiment if voyager is not installed
if (!requireNamespace("voyager", quietly = TRUE)) {
  if (!requireNamespace("SpatialFeatureExperiment", quietly = TRUE) ||
      packageVersion("SpatialFeatureExperiment") < "1.7.3") {
    message("Installing SpatialFeatureExperiment (>= 1.7.3)...")
    BiocManager::install("SpatialFeatureExperiment")
  } else {
    message("SpatialFeatureExperiment version is sufficient (>= 1.7.3), skipping installation.")
  }
} else {
  message("Voyager is installed. Skipping SpatialFeatureExperiment version check.")
}

#  Install Voyager if not already installed
if (!requireNamespace("voyager", quietly = TRUE)) {
  BiocManager::install("Voyager")
}

#  Install GitHub dependencies
install_github_if_missing("dmcable/spacexr", build_vignettes=FALSE)
install_github_if_missing("Center-for-Spatial-OMICs/SpatialQM")

#  Resolve conflict between dplyr and gridExtra
suppressWarnings({
  library(dplyr)
  library(gridExtra)
})

message("Dependency installation complete!")

