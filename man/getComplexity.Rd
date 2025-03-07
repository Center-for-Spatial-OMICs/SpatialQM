% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{getComplexity}
\alias{getComplexity}
\title{getComplexity}
\usage{
getComplexity(
  seu_obj = NULL,
  features = NULL,
  expMat = "path_to_expMat",
  platform = NULL
)
}
\arguments{
\item{seu_obj}{An optional Seurat object containing the expression data. If provided,
the function will calculate complexity based on the data in the object.}

\item{features}{Optional character vector of features to focus the complexity calculation
on a subset of all features. If NULL, all features are considered.}

\item{expMat}{The path to the expression matrix file or directory containing
expression data files. This parameter is used if no Seurat object is provided.
This function supports processing of '.mtx.gz' for sparse matrices, or '.h5'
files for HDF5 format expression matrices.}

\item{platform}{The platform from which the data originates. Valid options are 'Xenium', 'CosMx',
and 'Merscope'.}
}
\description{
Calculate the complexity of expression data, determining the number
of features (genes) required to reach half of the total expression sum in a dataset.
}
\details{
Calculate the total expression and determines the minimal set of features accounting for
at least half of this total, either across all features or a specified subset. If a Seurat
object is used, it derives the complexity from the RNA assay counts.
}
