# MinNetRank

## Install:

```r
install.packages("devtools")

library(devtools)

install_github("weitinging/MinNetRank")

install.packages("MinNetRank")

library(MinNetRank)
```

## example:

```r
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
```
