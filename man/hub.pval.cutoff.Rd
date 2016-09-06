\name{hub.pval.cutoff}
\alias{hub.pval.cutoff}
\title{Identification of hub genes in a Gene Co-expression Network based on gene connection significance values}
\description{The function returns list of top ranked hub genes under a particular condition/trait}
\usage{hub.pval.cutoff(x, beta, m, s, n)}
\arguments{
  \item{x}{x is a data frame of gene expression values where rows represent genes and columns represent samples/time point under a particular condition/trait.}
  \item{beta}{beta is a soft threshold parameter determined from the scale free property of biological networks (Gene Co-expression Network).}
  \item{m}{m is a scalar representing sample size and is less than or equal to number of columns in x.}
  \item{s}{s is a scalar representing number of times  each of the m samples will be resampled.}
  \item{n}{n is a scalar representing number of top ranked hub genes to be selected under that particular condition/trait}
}
\value{The function returns a vector containing list of top ranked hub genes for that particular condition/trait}
\author{Samarendra Das and Baidya Nath Mandal}
\examples{
data(rice_salt)
x=as.data.frame(rice_salt)
beta=6
m=18
s=80
n=20
hub.pval.cutoff(x, beta, m, s, n)
}
\keyword{gene}
\keyword{co-expression}
\keyword{network}
\keyword{hub gene}
\keyword{top rank}
\keyword{p value}