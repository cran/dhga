\name{hub.wgs}
\alias{hub.wgs}
\title{Identification of hub genes in a Gene Co-expression Network  based on Weighted Gene Score}
\description{The function returns the list of hub genes in a Gene Co-expression Network based on Weighted Gene Score without considering statistical significance values.}
\usage{hub.wgs(x, beta, n)}
\arguments{
  \item{x}{x is a data frame of gene expression values where rows represent genes and columns represent samples/time.}
  \item{beta}{beta is a soft threshold parameter determined from the scale free property of biological networks (Gene Co-expression Network).}
  \item{n}{n is a scalar representing number of genes to be declared as hub in the Gene Co-expression Network based on Weighted Gene Score.}
}
\value{The function returns a list of hub genes in the Gene Co-expression Network}
\author{Samarendra Das and Baidya Nath Mandal}
\examples{
data(rice_salt)

hub.wgs(rice_salt,beta=6, n=20)
}
\keyword{gene}
\keyword{co-expression}
\keyword{network}
\keyword{hub gene}
\keyword{Weighted Gene Score}