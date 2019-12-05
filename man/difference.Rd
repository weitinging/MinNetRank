\name{difference}
\alias{difference}
\title{
minimizing abs(restarted - diffused)
}
\description{
We choose β to minimize the absolute of the difference between diffusion and the retainment
}
\usage{
difference(beta, A_Matrix = AdjacencyMatrix)
}
\arguments{
  \item{beta}{the restart probability}
  \item{A_Matrix}{the adjacency matrix of interaction network}
}
\value{
%%  ~Describe the value returned
numeric
}
\author{
Ting Wei <weitinging@sjtu.edu.cn>; Zhangsheng Yu
Maintainer: Ting Wei <weitinging@sjtu.edu.cn>
}