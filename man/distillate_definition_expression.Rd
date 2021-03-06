% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distillate_definition_expression.R
\name{distillate_definition_expression}
\alias{distillate_definition_expression}
\title{Calculates the expression-fullness of each hierarchical level of a metabolic pathway/module}
\usage{
distillate_definition_expression(
  sample,
  definition_expression,
  def_table,
  level,
  expression_table
)
}
\arguments{
\item{definition_expression}{Definition-expression string}

\item{def_table}{Decomposed hierarchy matrix produced by create_step_matrix.R}

\item{level}{Hierarchical level}

\item{expression_vector}{Vector of expression values of functional units present in the genome}
}
\value{
A (partially) distilled definition string
}
\description{
Calculates the expression-fullness of each hierarchical level of a metabolic pathway/module
}
\examples{
distillate_definition_expression(definition_expression, def_table, level, present)
}
