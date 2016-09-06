\name{pvalue.hub}
\alias{pvalue.hub}
\title{Compuation of  gene connection sigificance values in a Gene Co-expression network}
\description{The function computes ths statisical significance values of gene connections in a GCN for a particular condition/trait/stress}
\usage{pvalue.hub(x, beta, m, s, plot=TRUE)}
\arguments{
  \item{x}{x is a data frame of gene expression values where rows represent genes and columns represent samples/time point under a particular condition/trait/stress.}
  \item{beta}{beta is a soft threshold parameter determined from the scale free property of biological networks (Gene Co-expression network).}
  \item{m}{m is a scalar representing sample size and is less than or equal to number of columns in x.}
  \item{s}{s is a scalar representing number of times  each of the m samples will be resampled.}
\item{plot}{plot is a character must either take logical value TRUE/FALSE representing whether the plot of the pvalues of gene connections of all genes in Gene Co-expression Network needs to be constructed or not.}
}
\value{The function returns a vector of p-values for all the genes in the GCN}
\author{Samarendra Das and Baidya Nath Mandal}
\examples{
data(rice_salt)
x=as.data.frame(rice_salt)
beta=6
m=18
s=80
pvalue.hub(x, beta, m, s, plot=TRUE)
}
\keyword{gene}
\keyword{co-expression}
\keyword{network}
\keyword{hub gene}
\keyword{p value}
\keyword{plot}