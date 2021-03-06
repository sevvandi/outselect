% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Feature_Functions.R
\name{ComputeMetaFeaturesMM}
\alias{ComputeMetaFeaturesMM}
\title{Computes features from a dataset, but using only MIN_MAX normalization.}
\usage{
ComputeMetaFeaturesMM(dat.o)
}
\arguments{
\item{dat.o}{This is the data with the last column named as \code{outlier} having values \code{yes} or \code{no}. All other columns are numeric. Non-numerical data need to be converted to numerical before computing features.}
}
\value{
The set of features.
}
\description{
Computes features from a dataset, but using only MIN_MAX normalization.
}
\examples{
data(Arrhythmia_withoutdupl_05_v05)
dat <- Arrhythmia_withoutdupl_05_v05
feat <- ComputeMetaFeaturesMM(dat)

}
