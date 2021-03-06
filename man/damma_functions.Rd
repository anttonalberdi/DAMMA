% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dama_functions.R
\name{damma_functions}
\alias{damma_functions}
\title{Aggregates compound values into indices at the function level}
\usage{
damma_functions(compounds_table, functions_table, transform = TRUE)
}
\arguments{
\item{transform}{Whether to transform to 0-1 scale. Default=TRUE}

\item{compound_table}{Compound table outputted by aggregate_compounds()}

\item{functions}{Table containing definitions and metadata of metabolic functions (provided by GAMMA)}
}
\value{
A function table aggregated at the compound level
}
\description{
Aggregates compound values into indices at the function level
}
\examples{
dama_functions(compound_table)
}
