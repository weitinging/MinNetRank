\name{MinNetRank}
\alias{MinNetRank}
\title{
MinNetRank
}
\description{
Discovery of candidate cancer genes by network-based integration of multi-omics data
}
\usage{
MinNetRank(Network = "AdjacencyMatrix", SNP = FALSE, TumorExpression = FALSE, NormalExpression = FALSE, CGC = KnownGenes, beta = 0.4841825)
}
\arguments{
  \item{Network}{the interaction network}
  \item{SNP}{the mutation matrix}
  \item{TumorExpression}{the tumor sample expresion}
  \item{NormalExpression}{the normal sample expresion}
  \item{CGC}{the known cancer genes}
  \item{beta}{the restart probability}
}
\details{
MinNetRank utilized minimum strategy to prioritize genes both the mutation relevance score and expression relevance score are high. MinNetRank was a single sample network diffusion approach that could detect personalized driver genes. MinNetRank combined the ranking of genes for individual samples into a robust population-level gene ranking. 
}
\value{
%%  ~Describe the value returned
%%  \item{Gene}{'the gene names'}
%%  \item{Score}{'the score of gene'}
%%  \item{KnownCancerGenes}{'this gene is known cancer gene (in CGC dataset) or not'}
%%  \item{Freq}{'the gene mutation frequency'}
}

\author{
Ting Wei <weitinging@sjtu.edu.cn>; Zhangsheng Yu
Maintainer: Ting Wei <weitinging@sjtu.edu.cn>
}

\examples{
library(MinNetRank)

#load the adjacency network
data("AdjacencyMatrix")

#load the known cancer genes
data("KnownGenes")

#load the mutation data
data("LihcMutation")

#load the tumor expression data
data("LihcTumorExpression")

#load the normal expression data
data("LihcNormalExpression")

#Using AdjacencyMatrix
##Using the mutation and expression data
Network = "AdjacencyMatrix"
LihcMinNetRank = MinNetRank(Network, SNP=LihcMutation, TumorExpression=LihcTumorExpression, NormalExpression=LihcNormalExpression, CGC=KnownGenes, beta = 0.4841825)
write.table(LihcMinNetRank, file='TCGA-LIHC.MinNetRank.Result.xls', quote =F, sep="\t", row.names = FALSE)

}