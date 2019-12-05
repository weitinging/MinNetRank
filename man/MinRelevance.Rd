\name{MinRelevance}
\alias{MinRelevance}
\title{
the minimum of MutationRelevanceScores and ExpressionRelevanceScores
}
\description{
A large score means that one gene is proximal to many mutation genes and differentially expressed genes
}
\usage{
MinRelevance(MutationRelevanceScores, ExpressionRelevanceScores)
}
\arguments{
  \item{MutationRelevanceScores}{the mutation relevance scores }
  \item{ExpressionRelevanceScores}{the expression relevance scores}
}
\value{
%%  ~Describe the value returned
matrix
}
\author{
Ting Wei <weitinging@sjtu.edu.cn>; Zhangsheng Yu
Maintainer: Ting Wei <weitinging@sjtu.edu.cn>
}
