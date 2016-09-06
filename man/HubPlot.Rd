\name{HubPlot}
\alias{HubPlot}
\title{Venn plot of hub genes under two conditions namely stress and control}
\description{The function produces a venn diagram of number of hub genes and unique hub genes under two conditions}
\usage{HubPlot(pvalue.stress, pvalue.control, alpha)}
\arguments{
  \item{pvalue.stress}{pvalue.stress is a vector of p-values for genes in the GCN under stress condition. pvalue.stress can be obtained from the output of hub.stress function.}
  \item{pvalue.control}{pvalue.control is a vector of p-values for genes in the GCN under normal or control condition. pvalue.control can be obtained from the output of hub.control function.}
  \item{alpha}{alpha is a scalar representing statistical level of significance. Default is alpha=0.0001}
}
\value{The function produces a venn plot of hub genes and unique hub genes under two conditions}
\author{Samarendra Das and Baidya Nath Mandal}
\examples{
data(rice_salt)
pval.stres <- pvalue.hub(rice_salt, beta=6, m=18, s=80, plot=FALSE)
pvalue.stress <- pval.stres[, 2]

data(rice_normal)
pval.control <- pvalue.hub(rice_normal, beta=6, m=18, s=80, plot=FALSE)
pvalue.control <- pval.control[, 2]
HubPlot(pvalue.stress, pvalue.control, alpha=0.0001)
}
\keyword{gene}
\keyword{co-expression}
\keyword{network}
\keyword{hub gene}
\keyword{venn plot}
\keyword{level of significance}