\name{WeightedGeneScore}
\alias{WeightedGeneScore}
\title{Computation of weighted gene score of genes in Gene Co-expression Network}
\description{The function computes the weighted gene scores of all genes in a Gene Co-expression Network}
\usage{WeightedGeneScore(x, beta, plot=TRUE)}
\arguments{
  \item{x}{x is a data frame of gene expression values where rows represent genes and columns represent samples/time points under a particular condition/stress/trait.}
  \item{beta}{beta is a soft threshold parameter determined from the scale free property of biological networks (Gene Co-expression Network).}
  \item{plot}{plot is a character must either take logical value TRUE/FALSE representing whether the plot of the gene connection degrees of all genes in Gene Co-expression Network needs to be constructed or not.}
}
\value{The function returns a vector of weighted gene scores of genes (nodes) in the Gene Co-expression Network along with the gene connectivity plot.}
\author{Samarendra Das and Baidya Nath Mandal}
\examples{
data(rice_normal)
beta=6
WeightedGeneScore(rice_normal, beta, plot=FALSE)
}
\keyword{gene}
\keyword{co-expression}
\keyword{network}
\keyword{weighted gene score}
\keyword{plot}