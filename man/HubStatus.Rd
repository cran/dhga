\name{HubStatus}
\alias{HubStatus}
\title{Hub status of the genes in a gene co-expression network}
\description{The function returns hub status of each gene in a gene co-expression network (GCN)}
\usage{HubStatus(pvalue.stress, pvalue.control, alpha)}
\arguments{
  \item{pvalue.stress}{pvalue.stress is a vector of p-values for genes in the GCN under stress condition. pvalue.stress can be obtained from the output of hub.stress function.}
  \item{pvalue.control}{pvalue.control is a vector of p-values for genes in the GCN under normal or control condition. pvalue.control can be obtained from the output of hub.control function.}
  \item{alpha}{alpha is a scalar representing statistical level of significance. Default is alpha=0.0001}
}
\value{The function returns a list with two components. First component returns the list of genes along with their hub status. Second component gives a table containing number of genes under different categories of hubs viz. housekeeping, unique to stress and unique to normal}
\author{Samarendra Das and Baidya Nath Mandal}
\examples{
data(rice_salt)
pval.stres <- pvalue.hub(rice_salt, beta=6, m=18, s=80, plot=FALSE)
p1 <- pval.stres[, 2]

data(rice_normal)
pval.control <- pvalue.hub(rice_normal, beta=6, m=18, s=80, plot=FALSE)
p2 <- pval.control[, 2]
HubStatus(p1,p2,alpha=0.0001)
}
\keyword{gene}
\keyword{co-expression}
\keyword{network}
\keyword{hub gene}
\keyword{hub status}
\keyword{level of significance}