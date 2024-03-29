% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_MCI_community.R
\name{compute_MCI_community}
\alias{compute_MCI_community}
\title{Calculates the MCI of a metabolic pathway at the community level}
\usage{
compute_MCI_community(definition, abundance_vector)
}
\arguments{
\item{definition}{Definition string of a given metabolic pathway/module}

\item{abundance_vector}{Vector containing relative abundance data per identifier}
}
\value{
A MCI value
}
\description{
Calculates the MCI of a metabolic pathway at the community level
}
\examples{
compute_MCI_community(definition,abundance_vector)
}
