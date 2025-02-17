% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{getAllMetrics}
\alias{getAllMetrics}
\title{getAllMetrics}
\usage{
getAllMetrics(df_samples, features = NULL)
}
\arguments{
\item{df_samples}{A data.frame that must contain the following columns:
\itemize{
\item{sample_id}{Unique identifier for each sample.}
\item{platform}{The platform used, either "Xenium" or "CosMx".}
\item{expMat}{Path to the expression matrix file. For Xenium platform, this should be a path to the cell_feature_matrix folder.}
\item{tx_file}{Path to the transcript file.}
\item{cell_meta}{Path to the cell metadata file.}
}}
}
\value{
A modified version of df_samples with the computed metrics added as new columns.
}
\description{
Compute various metrics for a given set of samples
}
\details{
Computes various metrics for spatial transcriptomic data for a given set of samples. This function processes a data.frame containing sample information, including paths to necessary files, to compute a variety of metrics for each sample.
The metrics outputted include the number of cells, specificity FDR, number of transcripts per cell, transcripts per area of segmented cell, transcripts per nucleus, signal to noise ratio, fraction of transcripts in segmented cells, Mutually Exclusive Co-Expression Rate (MECR), sparsity, and Shannon entropy.
To initiate vignette, the data frame must include the columns: sample_id, platform, expMat, tx_file, and cell_meta.
}
