% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/damma.R
\name{damma}
\alias{damma}
\title{Generate the fullness table of metabolic pathways/modules from a MAG annotation table}
\usage{
damma(annotations, functions, genomecol, keggcol, eccol, pepcol)
}
\arguments{
\item{annotations}{Table containing MAG identifiers and annotation codes}

\item{functions}{Table containing definitions and metadata of metabolic functions (provided by GAMMA)}

\item{genomecol}{Column index (number) of the annotations table containing the MAG identifiers}

\item{keggcol}{Column index(es) of the annotations table in which to search for KEGG KO annotations}

\item{eccol}{Column index(es) of the annotations table in which to search for Enzyme Commision (EC) annotations}

\item{pepcol}{Column index(es) of the annotations table in which to search for Peptidase annotations}
}
\value{
A fullness matrix
}
\description{
Generate the fullness table of metabolic pathways/modules from a MAG annotation table
}
\examples{
damma(annotations,functions,genomecol,keggcol,eccol,pepcol)
damma(annotations,functions,genomecol=2,keggcol=9,eccol=c(10,19),pepcol=12)
}
