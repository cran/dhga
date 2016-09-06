##############Dependent Packages####################################
library(VennDiagram)                                               #
library(stats)
library(graphics)
######################Ends Here#####################################

####################################################################
#                     Differential Hub Gene Analysis                #                                
####################################################################
Adjacency <- function (x, beta, threshold)
  
{
  this.call = match.call()
  if ((!class(x)=="data.frame")) 
  {
    warning("x must be a data frame and rows as gene names")
  }
  if (!(class(beta)=="numeric" & beta > 1 )) 
  {
    stop ("beta must be integer and  must be beta > 1")
  } 
  
  if (missing (beta))
  {
    beta <- 6
  }
  if (!(class(threshold)=="numeric")) 
  {
    stop ("threshold must be numeric")
  } 
  if (threshold <= 0 & threshold >= 1) 
  {
    stop ("threshold must be any value between 0 and 1")
  } 
  gene <- rownames(x)
  x <- as.matrix(x)
  M <- ncol(x)    #####number of samples
  n <- nrow(x)    #####number of genes
  Adjacency <- ((abs(cor(t(x), method="pearson", use="p")))^beta)-diag(n)
  rownames(Adjacency) <- gene
  colnames(Adjacency) <- gene
  
  Adj1 <- ifelse(Adjacency >= threshold, 1, 0)
  Adj1[lower.tri(Adj1)] <- 0
  idx <- which(Adj1 == 1)
  idx.m1 <- idx -1
  rows <- idx.m1 %% nrow (Adj1) + 1
  cols <- idx.m1 %/% nrow (Adj1) + 1
  EdgeList <- data.frame(Node_i=rownames(Adj1)[rows],
                                 Node_j=colnames(Adj1)[cols],
                                 Source=rep('Predicted Interaction', length(rows)),
                                 stringsAsFactors=FALSE)
  NodeList <- union(rownames(Adj1)[rows], colnames(Adj1)[cols])
  NodeList <- as.data.frame(NodeList)
  out <- list(Adjacency=Adjacency, NodeList=NodeList, EdgeList=EdgeList)
  class(out) <- c("Adjacency matrix", "genelist", "edgelist")
  return (out)
}

#########################Weighted gene score###################################
WeightedGeneScore <- function(x, beta, plot=TRUE)
  {
     this.call = match.call()
     if ((!class(x)=="data.frame")) 
     {
       warning("x must be a data frame and rows as gene names")
     }
     if (!(class(beta)=="numeric" & beta > 1 )) 
     {
       stop ("beta must be integer and  must be beta > 1")
     } 
     gene <- rownames (x)
     x <- as.matrix(x)
     M <- ncol(x)    #####number of samples
     n <- nrow(x)    #####number of genes
     Adj.pop <- ((abs(cor(t(x), method="pearson", use="p")))^beta)-diag(n)
      if (missing (beta))
        {
        beta <- 6
        Adj.pop <- ((abs(cor(t(x), method="pearson", use="p")))^beta)-diag(n)
        }
        WGS <- as.vector(apply(Adj.pop, 2, sum, na.rm=T))
        WGS <- as.data.frame(WGS)
        rownames(WGS) <- gene
      if (plot==TRUE)
      {
        wgs <- as.vector(WGS[,1])
        ranks <- sort(wgs, decreasing=TRUE, index.return=TRUE)$ix
        wgs_sort <- wgs[ranks]
        plot(1:length(ranks), wgs_sort, main="Weighted Gene Score Fig", type="p", xlab="Genes", ylab="WGS")
      }
     #class(WGS) <- "weighted gene connectivity"
     WGS
}
     
 #################hub gene identification based on WGS##############    
hub.wgs <- function (x, beta, n)
  
{
  this.call = match.call()
  if ((!class(x)=="data.frame")) 
  {
    warning("x must be a data frame and rows as gene names")
  }
  if (!(class(beta)=="numeric" & beta > 1 )) 
  {
    stop ("beta must be integer and  must be beta > 1")
  } 
  
  if (missing (beta))
  {
    beta <- 6
  }
  if(!(class(n) =="numeric" & n < nrow(x)))
  {
    stop("n must be numeric and less than the total number of genes")
  }
  gene <- rownames(x)
  score <- as.vector(WeightedGeneScore(x, beta, plot=FALSE)$WGS)
  id.cond <- sort(-score, index.return=TRUE)$ix
  id.cond1 <- id.cond[1:n]
  select.hub.cond1 <- gene[id.cond1]
  class(select.hub.cond1) <- "Hub Genes"
  return (select.hub.cond1)
}

#####################Hub identification based on statistical significance value under stress###################

pvalue.hub <- function (x, beta, m, s, plot=TRUE)
{  ####x must be data frame####
   this.call = match.call()
   if ((!class(x)=="data.frame")) 
   {
     warning("x must be a data frame and rows as gene names")
   }
   if (!(class(beta)=="numeric" & beta > 1 )) 
   {
     stop ("beta must be integer and  must be beta > 1")
   } 
   if (!(class(m)=="numeric" & m <= ncol(x))) 
   {
     stop ("m must be integer and less than total number of samples")
   }
   if (!(class(s)=="numeric" & s > 30)) 
   {
     warning("s must be integer and sufficiently large")
   }
   if (missing (beta))
   {
     beta <- 6
   }
   gene <- rownames (x)
   x <- as.matrix(x)
   M <- ncol(x)    #####number of samples
   n <- nrow(x)    #####number of genes
   Adj.pop <- ((abs(cor(t(x), method="pearson", use="p")))^beta)-diag(n)
   conn.tot <- as.vector(apply(Adj.pop, 2, sum, na.rm=T))
   pop.mean <- as.numeric(mean (conn.tot))  ####Average connection degree of the complete network model
   #############Resampling Procedure ###################################
   Deg <- matrix(0, n, s)
   for (i in 1: s){
     samp <- sample(M, m, replace=T)
     dat.samp <- x[, samp]
     ADJ.samp <- ((abs(cor(t(dat.samp), method="pearson", use="p")))^beta)-diag(n)
     con.samp <- as.vector(apply(ADJ.samp, 2, sum, na.rm=T))
     Deg[,i] <- con.samp
   }
   mu <- round(pop.mean, 2)                       # value under null hypothesis
   D <- Deg - mu                                  # Transformed score of adjacency scores under H0
   pval.deg <- vector(mode="numeric", length=n)
   for (i in 1: n) {
     Z <- D[i,]
     Z <- Z[Z != 0]
     n1 <- length(Z)
     r <- rank(abs(Z))
     tplus <- sum(r[Z > 0])
     etplus <- n1 * (n1 + 1) / 4
     vtplus <- n1 * (n1 + 1) * (2 * n1 + 1) / 24
     p.value=pnorm(tplus, etplus, sqrt(vtplus), lower.tail=FALSE)
     pval.deg[i]=p.value
   }
   res <- cbind(conn.tot, pval.deg)
   rownames (res) <- gene
   colnames (res) <- c("Total degree", "P value")
   if (plot==TRUE)
   {
     pval <- as.vector(res[,2])
     ranks <- sort(pval, decreasing=FALSE, index.return=TRUE)$ix
     pval_sort <- pval[ranks]
     plot(1:length(ranks), pval_sort, main="Hub gene selection plot", type="p", xlab="Genes", ylab="pvalue for WGS")
   }
   #out1 <- list(conn.tot=conn.tot, pval.deg=pval.deg)
   class(res) <- c("Statistical significance of gene connections")
   res
}


########################Identification of list of hubs under stress###################

hub.pval.cutoff <- function (x, beta, m, s, n)
{
  this.call = match.call()
  if ((!class(x)=="data.frame")) 
  {
    warning("x must be a data frame and rows as gene names")
  }
  if (!(class(beta)=="numeric" & beta > 1 )) 
  {
    stop ("beta must be integer and  must be beta > 1")
  } 
  if (!(class(m)=="numeric" & m <= ncol(x))) 
  {
    stop ("m must be integer and less than total number of samples")
  }
  if (!(class(s)=="numeric" & s > 30)) 
  {
    warning("s must be integer and sufficiently large")
  }
  if (missing (beta))
  {
    beta <- 6
  }
  if(!(class(n) =="numeric" & n < nrow(x)))
  {
    stop("n must be numeric and less than the total number of genes")
  }
  gene <- rownames(x)
  pval <- pvalue.hub(x, beta, m, s, plot=FALSE)
  pvalue <- as.vector(pval[, 2])
  id.cond <- sort(pvalue, decreasing= FALSE, index.return=TRUE)$ix
  id.cond1 <- id.cond[1:n]
  select.hub.cond1 <- gene[id.cond1]
  class(select.hub.cond1) <- "Hub Genes"
  return (select.hub.cond1)
}


########Venn plotting############

HubPlot <- function (pvalue.stress, pvalue.control, alpha)
{
  if(!(class(alpha)=="numeric" & alpha < 1))
  {
    stop ("Alpha (Level of significance) must be less than 1")
  }
  id1 <- which (pvalue.stress < alpha)
  id2 <- which (pvalue.control < alpha)
  hub1 <- names (pvalue.stress[id1])
  hub2 <- names (pvalue.control[id2])
  hub.common <- intersect (hub1, hub2)
  venn.plot <- draw.pairwise.venn(area1=length(hub1), 
                                  area2=length(hub2), 
                                  cross.area=length(hub.common),
                                  category = c("Stress", "Control"),
                                  fill = c("red", "green"), lty = "solid", cex = 2, 
                                  cat.cex = 1,cat.col=c("red", "green"))
  class (venn.plot) <- "Hub Plot"
  #return (venn.plot)                        
}

#########################Hub Status####################

HubStatus <- function(pvalue.stress, pvalue.control, alpha)
{
  if(!(class(alpha)=="numeric" & alpha < 1))
  {
    stop ("Alpha (Level of significance) must be less than 1")
  }
  genes <- names(pvalue.stress)
  status <- ifelse(pvalue.stress < alpha & pvalue.control < alpha,"Housekeeping Hub",
                   ifelse(pvalue.stress < alpha & pvalue.control > alpha, "Unique Hub to Stress",
                          ifelse(pvalue.stress > alpha & pvalue.control < alpha,"Unique Hub to Normal","Not Hub")))
  hub.status <- as.data.frame(cbind(genes, status), row.names=FALSE, sep="\t")
  colnames(hub.status) <- c("Genes", "Hub Status")
  out1 <- table(hub.status[,2])
  
  out11 <- list(hub.status=hub.status, out1=out1)
  class(out11) <- c("HubStatus_Genes", "HubStatus_Matrix")
  out11
}     

#######################Differential Hub Analysis ###############################

DiffHub <-function(x,y,m1,m2,s,beta,alpha, plot=TRUE)
{
  this.call = match.call()
  if (!(class(x)=="data.frame" & class(y)=="data.frame")) 
  {
    warning("x must be a data frame and rows as gene names")
  }
  if (!(class(beta)=="numeric" & beta > 1 )) 
  {
    stop ("beta must be integer and  must be beta > 1")
  } 
  if (!(class(m1)=="numeric" & m1 <= ncol(x))) 
  {
    stop ("m1 must be integer and less than total number of samples")
  }
  if (!(class(m2)=="numeric" & m2 <= ncol(y))) 
  {
    stop ("m2 must be integer and less than total number of samples")
  }
  if (!(class(s)=="numeric" & s > 30)) 
  {
    warning("s must be integer and sufficiently large")
  }
  if (missing (beta))
  {
    beta <- 6
  }
  if(!(class(alpha)=="numeric" & alpha < 1))
  {
    stop ("Alpha (Level of significance) must be less than 1")
  }
  temp <- match(rownames(x), rownames(y))
  if (sum(is.na(temp)) == length(temp)) {
    stop("The genes in x and y must match with each other")
  }
  genes <- rownames(x)
  pval.stres <- pvalue.hub(x, beta, m1, s, plot=FALSE)
  pvalue.stress <- pval.stres[, 2]
  pval.control <- pvalue.hub(y, beta, m2, s, plot=FALSE)
  pvalue.control <- pval.control[, 2]
  if (plot=="TRUE")
  {
    HubPlot(pvalue.stress,pvalue.control,alpha)
  }
  pvalue.stress1 <- as.vector(pval.stres[, 2])
  pvalue.control1 <- as.vector(pval.control[, 2])
  status <- ifelse(pvalue.stress1 < alpha & pvalue.control1 < alpha,"Housekeeping Hub",
                   ifelse(pvalue.stress1 < alpha & pvalue.control1 > alpha, "Unique Hub to Stress",
                          ifelse(pvalue.stress1 > alpha & pvalue.control1 < alpha,"Unique Hub to Normal","Not Hub")))
  out <- as.data.frame(cbind(genes, pvalue.stress1,pvalue.control1, status), row.names=FALSE, sep="\t")
  colnames(out) <- c("Genes", "pval_stres", "pval_control", "HubStatus")
  #class(out) <- "Differental Hub Analysis"
  out
}  

#########################Ends here####################
#Execute
#load("C:/My package/dhga/data/rice_salt.rda")
#load("C:/My package/dhga/data/rice_normal.rda") 
#DiffHub(dat.trt,dat.ctrl,10,12,50,beta=6,alpha=0.0001)
