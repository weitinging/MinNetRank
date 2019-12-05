\name{MinNetRank-package}
\alias{MinNetRank-package}
\alias{MinNetRank}
\docType{package}
\title{MinNetRank}
\description{
Discovery of candidate cancer genes by network-based integration of multi-omics data
}
\details{
MinNetRank utilized minimum strategy to prioritize genes both the mutation relevance score and expression relevance score are high. MinNetRank was a single sample network diffusion approach that could detect personalized driver genes. MinNetRank combined the ranking of genes for individual samples into a robust population-level gene ranking. 
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

##Using the mutation data
Network = "AdjacencyMatrix"
LihcMinNetRankMutation = MinNetRank(Network, SNP=LihcMutation, CGC=KnownGenes, beta = 0.4841825)
write.table(LihcMinNetRankMutation, file='TCGA-LIHC.MinNetRank.Mutation.xls', quote =F, sep="\t", row.names = FALSE)

##Using the expression data
Network = "AdjacencyMatrix"
LihcMinNetRankExp = MinNetRank(Network, TumorExpression=LihcTumorExpression, NormalExpression=LihcNormalExpression, CGC=KnownGenes, beta = 0.4841825)
write.table(LihcMinNetRankExp, file='TCGA-LIHC.MinNetRank.Expression.xls', quote =F, sep="\t", row.names = FALSE)
}