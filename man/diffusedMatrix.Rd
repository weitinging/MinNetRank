\name{diffusedMatrix}
\alias{diffusedMatrix}
\title{
diffusing the normalized adjacency matrix
}
\description{
We consider the direct and indirectly neighbors of each gene in the interaction network. The value reflects the connectivity between nodes j and i. 
}
\usage{
diffusedMatrix(M_Matrix, beta)
}
\arguments{
  \item{M_Matrix}{the normalized adjacency matrix}
  \item{beta}{the restart probability}
}
\value{
%%  ~Describe the value returned
matrix
}
\author{
Ting Wei <weitinging@sjtu.edu.cn>; Zhangsheng Yu
Maintainer: Ting Wei <weitinging@sjtu.edu.cn>
}
