% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{plotARI}
\alias{plotARI}
\title{plotARI}
\usage{
plotARI(cluster_metrics)
}
\arguments{
\item{cluster_metrics}{Data frame with \code{sample_id}, \code{value}, and \code{metric} columns, filtered for "ARI".}
}
\value{
A ggplot object depicting ARI values for each sample.
}
\description{
Generates a scatter plot visualizing Adjusted Rand Index (ARI) values across samples to assess clustering performance.
}
