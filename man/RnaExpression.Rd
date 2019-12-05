\name{RnaExpression}
\alias{RnaExpression}
\title{
geting the differential expression score for each individual sample
}
\description{
We canculate the normalized matrix of differential expression by column.
}
\usage{
RnaExpression(Expression, Network)
}
\arguments{
  \item{Expression}{differential expression with the paired tumor and normal expression}
  \item{Network}{the interaction network}
}
\value{
%%  ~Describe the value returned
matrix
}
\author{
Ting Wei <weitinging@sjtu.edu.cn>; Zhangsheng Yu
Maintainer: Ting Wei <weitinging@sjtu.edu.cn>
}