\name{Adjacency}
\alias{Adjacency}
\title{Construction of Gene Co-expression Network}
\description{The function computes the adjacency matrix using the soft threshold parameter to construct a Gene Co-expression Network}
\usage{Adjacency(x, beta, threshold)}
\arguments{
  \item{x}{x is a data frame of gene expression values where rows represent genes and columns represent samples/time points under a particular condition/stress/trait.}
  \item{beta}{beta is a soft threshold parameter determined from the scale free property of biological networks (Gene Co-expression Network).}
  \item{threshold}{threshold is a scalar representing the threshold value for the adjacency matrix to construct the Gene Co-expression Network.}
}
\value{The function returns a list of outputs like adjacency matrix, Node-list and Edge-list for constructing a Gene Co-expression Network at certain threshold level}
\author{Samarendra Das and Baidya Nath Mandal}
\examples{
data(rice_normal)
beta=6
threshold=0.4
Adjacency(rice_normal,6,0.4)
}
\keyword{gene}
\keyword{co-expression}
\keyword{network}
\keyword{Adjacency}
\keyword{Nodelist}
\keyword{Edgelist}