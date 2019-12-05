\name{Mutation}
\alias{Mutation}
\title{
geting the gene mutation score for each individual sample
}
\description{
SNP is the mutation matrix. If sample k has the gene i mutation, the value =1. We canculate the normalized matrix of SNP by column.
}
\usage{
Mutation(SNP, Network)
}
\arguments{
  \item{SNP}{the mutation matrix}
  \item{Network}{the interaction network}
}
\value{
%%  ~Describe the value returned
the normalized matrix of SNP by column.
}
\author{
Ting Wei <weitinging@sjtu.edu.cn>; Zhangsheng Yu
Maintainer: Ting Wei <weitinging@sjtu.edu.cn>
}

