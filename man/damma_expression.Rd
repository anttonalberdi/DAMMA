% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/damma_expression.R
\name{damma_expression}
\alias{damma_expression}
\title{Generate the fullness table of metabolic pathways/modules from a Genome annotation table}
\usage{
damma_expression(
  expression,
  annotations,
  functions,
  genecol,
  genomecol,
  keggcol,
  eccol,
  pepcol
)
}
\arguments{
\item{expression}{Table containing normalised gene expression data with genes in rows and samples in columns}

\item{annotations}{Table containing gene and genome identifiers, and annotation codes}

\item{functions}{Table containing definitions and metadata of metabolic functions (provided by GAMMA)}

\item{genecol}{Column index (number) of the annotations table containing the gene identifiers}

\item{genomecol}{Column index (number) of the annotations table containing the genome identifiers}

\item{keggcol}{Column index(es) of the annotations table in which to search for KEGG KO annotations}

\item{eccol}{Column index(es) of the annotations table in which to search for Enzyme Commision (EC) annotations}

\item{pepcol}{Column index(es) of the annotations table in which to search for Peptidase annotations}
}
\value{
A list of pathway-expression matrices (one table per genome)
}
\description{
Generate the fullness table of metabolic pathways/modules from a Genome annotation table
}
\examples{
damma_expression(expression,annotations,functions,genecol,genomecol,keggcol,eccol,pepcol)
damma_expression(expression,annotations,functions,genecol=1,genomecol=2,keggcol=9,eccol=c(10,19),pepcol=12)
}
