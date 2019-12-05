\name{Relevance}
\alias{Relevance}
\title{
the mutation relevance scores or the expression relevance scores
}
\description{
the product of D_Matrix and NormalizedScore
the value reflects the proximity of one gene to mutation genes or differential expression gene
}
\usage{
Relevance(D_Matrix, NormalizedScore)
}
\arguments{
  \item{D_Matrix}{the diffused matrix}
  \item{NormalizedScore}{the gene mutation score or the gene expression score for each individual sample}
}
\value{
%%  ~Describe the value returned
matrix
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Ting Wei <weitinging@sjtu.edu.cn>; Zhangsheng Yu
Maintainer: Ting Wei <weitinging@sjtu.edu.cn>
}
