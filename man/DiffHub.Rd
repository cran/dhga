\name{DiffHub}
\alias{DiffHub}
\title{Differential Hub status of the genes in a gene co-expression network}
\description{The function returns hub status of each gene in a gene co-expression network under two contrasting conditions}
\usage{DiffHub(x,y,m1,m2,s,beta,alpha, plot=TRUE)}
\arguments{
  \item{x}{x is a data frame of gene expression values where rows represent genes and columns represent samples/time point under stress condition.}
  \item{y}{y is a data frame of gene expression values where rows represent genes and columns represent samples/time points under control condition.}
  \item{m1}{m1 is a scalar representing sample size and is less than or equal to number of columns in x.}
  \item{m2}{m2 is a scalar representing sample size and is less than or equal to number of columns in y.}
  \item{s}{s is a scalar representing number of times  each of the m samples will be resampled.}
  \item{beta}{beta is a soft threshold parameter determined from the scale free property of biological networks (GCN).}
  \item{alpha}{alpha is a scalar representing statistical level of significance. Default is alpha=0.0001}
  \item{plot}{plot is a character representing whether the hubplot can be drawn (TRUE) or not (FALSE).}
}
\value{The function returns a list with two components. First component returns the list of genes along with their hub status. Second component gives a table containing number of genes under different categories of hubs viz. housekeeping, unique to stress and unique to normal. The function also returns a venn plot of hub genes and unique hub genes under two conditions.}
\author{Samarendra Das and Baidya Nath Mandal}
\examples{
data(rice_salt)
data(rice_normal)
DiffHub(rice_salt,rice_normal,18,18,80, 6, alpha=0.0001, plot=TRUE)
}
\keyword{gene}
\keyword{co-expression}
\keyword{network}
\keyword{hub gene}
\keyword{hub status}
\keyword{level of significance}