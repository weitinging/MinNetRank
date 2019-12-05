#1. Interaction network
##adjacency matrix A and the normalized adjacency matrix M
normalizedMatrix = function(A_Matrix){
  colsums = apply(A_Matrix, 2, sum)
  M_Matrix = t(t(A_Matrix)/(colsums+1e-16))
  return(M_Matrix)
}

##the diffused matrix D
diffusedMatrix = function(M_Matrix, beta){
  D_Matrix = beta*solve(diag(1, dim(M_Matrix)[1], dim(M_Matrix)[2]) - (1-beta)*M_Matrix)
  return(D_Matrix)
}

##we choose β to minimize |sum(D_ij [M_ij≠0])-sum(D_ij [M_ij=0])|
difference = function(beta, A_Matrix = AdjacencyMatrix){
  M_Matrix = normalizedMatrix(A_Matrix)
  D_Matrix = diffusedMatrix(M_Matrix, beta)
  restarted = sum(D_Matrix[A_Matrix != 0])
  diffused = sum(D_Matrix[A_Matrix = 0])  
  return(abs(restarted - diffused))
}

chooseBeta = function(){
  return(optimize(difference, interval = c(0.1, 0.9))$minimum)
} 

#2. Mutation data
##get the mutation score of each gene in each individual sample
Mutation = function(SNP, Network){
  SNPScore = matrix(0, dim(Network)[1], dim(SNP)[2])
  rownames(SNPScore) = rownames(Network)
  colnames(SNPScore) = colnames(SNP)
  SNPScore[intersect(rownames(SNPScore),rownames(SNP)), ] = as.matrix(SNP[intersect(rownames(SNPScore),rownames(SNP)),])
  GeneZero = colnames(SNPScore)[apply(SNPScore, 2, sum) == 0]
  SNPScore = SNPScore[, setdiff(colnames(SNPScore), GeneZero)]
  SNPScore = abs(SNPScore)
  NormalizedMutationScore = t(t(SNPScore)/apply(SNPScore, 2, sum))
  return(NormalizedMutationScore)
}

#3. Expresion data
##calculat the Absolute of Log2 Fold-Change (ALFC) of gene expression between the paired tumor and normal samples
TumorNormalExp = function(TumorExpression, NormalExpression ){
  ##get gene expression for genes if the tumor with paired normal
  matchedSample = intersect(colnames(TumorExpression),colnames(NormalExpression))
  matchedDiffExpression = log2(TumorExpression[,matchedSample] / NormalExpression[,matchedSample])
  
  #get the gene expression for genes if the tumor with no paried normal
  unmatchedSample = setdiff(colnames(TumorExpression) , matchedSample)
  generalExpression = apply(NormalExpression,1,mean)
  unmatchedDiffExpression = log2(TumorExpression[,unmatchedSample] / generalExpression)
  
  #combine the  expressions
  Expression<-cbind(matchedDiffExpression,unmatchedDiffExpression)
  Expression[is.na(Expression)] = 0
  Expression[Expression == "Inf"] = 0
  Expression[Expression == "-Inf"] = 0
  return(Expression)
}

##get the expression score of each gene in each individual sample
RnaExpression = function(Expression, Network){
  ExpressionScore = matrix(0, dim(Network)[1], dim(Expression)[2])
  rownames(ExpressionScore) = rownames(Network)
  colnames(ExpressionScore) = colnames(Expression)
  ExpressionScore[intersect(rownames(ExpressionScore),rownames(Expression)), ] = as.matrix(Expression[intersect(rownames(ExpressionScore),rownames(Expression)),])
  ExpressionScore = abs(ExpressionScore)
  NormalizedExpressionScore = t(t(ExpressionScore)/apply(ExpressionScore, 2, sum))
  return(NormalizedExpressionScore)
}

#4. the relevance  scores
Relevance = function(D_Matrix, NormalizedScore){
  ##the relevance  scores of mutation data or expression data
  RelevanceScores = D_Matrix%*%NormalizedScore
  return(RelevanceScores)
}

#5. the pmin of expression relevance  scores and mutation relevance  score
MinRelevance = function(MutationRelevanceScores, ExpressionRelevanceScores){
  ##the relevance  score of gene is 0 in all samples
  for (i in 1:dim(MutationRelevanceScores)[1]){
    if(sum(MutationRelevanceScores[i, ]==0)/dim(MutationRelevanceScores)[2] >0.95){
      MutationRelevanceScores[i, ] = 0
    }
  }
  for (i in 1:dim(ExpressionRelevanceScores)[1]){
    if(sum(ExpressionRelevanceScores[i, ]==0)/dim(ExpressionRelevanceScores)[2] >0.95){
      ExpressionRelevanceScores[i, ] = 0
    }
  }
  GeneZero = union(rownames(MutationRelevanceScores[(apply(MutationRelevanceScores,1,sum) == 0),]), rownames(ExpressionRelevanceScores[(apply(ExpressionRelevanceScores,1,sum) == 0),]))
  ExpressionRelevanceScores = ExpressionRelevanceScores[setdiff(rownames(MutationRelevanceScores), GeneZero), ]
  MutationRelevanceScores = MutationRelevanceScores[setdiff(rownames(MutationRelevanceScores), GeneZero), ]
  
  ##final relevance  scores
  samesamples = intersect(colnames(MutationRelevanceScores),colnames(ExpressionRelevanceScores))
  MinRelevanceScores = pmin(ExpressionRelevanceScores[ ,samesamples], MutationRelevanceScores[ ,samesamples])
  return(MinRelevanceScores)
}


#6.aggregate gene ranking in each individual sample to a robust population-level gene ranking
RankSum = function(RelevanceScores){
  geneRank = matrix(0, dim(RelevanceScores)[1], dim(RelevanceScores)[2])
  for (i in 1:dim(RelevanceScores)[2]){
    geneRank[, i] = (dim(RelevanceScores)[1] + 1) - rank(RelevanceScores[, i])
  }
  rownames(geneRank) = rownames(RelevanceScores)
  colnames(geneRank) = colnames(RelevanceScores)
  RankResult = rowSums(geneRank)
  RankResult = sort(RankResult)
  return(RankResult)
}

#7.Result
MinNetRank = function(Network = "AdjacencyMatrix", SNP = FALSE, TumorExpression  = FALSE, NormalExpression  = FALSE, CGC = KnownGenes, beta = 0.4841825){
  if (Network == "DiffusedMatrix"){
    Network = DiffusedMatrix
    print("Using DiffusedMatrix in MinNetWork.")
    if (identical(SNP, FALSE) & identical(TumorExpression, FALSE)){
      print("You should provide Mutation data or Expression data.")
    }
    else if(identical(SNP, FALSE)){
      print("Mutation data are not provided, Only using Expression data")
      RankResult = RankSum(Relevance(Network, RnaExpression(TumorNormalExp(TumorExpression, NormalExpression), Network)))
    }
    else if(identical(TumorExpression, FALSE)){
      print("Expression data are not provided, Only using Mutation data")
      RankResult = RankSum(Relevance(Network, Mutation(SNP, Network)))
    }
    else{
      print("Using Mutation data and Expression data.")
      RankResult = RankSum(MinRelevance(Relevance(Network, Mutation(SNP, Network)), Relevance(Network, RnaExpression(TumorNormalExp(TumorExpression, NormalExpression), Network))))  
    }
  }  
  else if (Network == "AdjacencyMatrix"){
    Network = AdjacencyMatrix
    print("Using AdjacencyMatrix in MinNetWork.")
    if (identical(SNP, FALSE) & identical(TumorExpression, FALSE)){
      print("You should provide Mutation data or Expression data.")
    }
    else if(identical(SNP, FALSE)){
      print("Mutation data are not provided, Only using Expression data")
      RankResult = RankSum(Relevance(diffusedMatrix(normalizedMatrix(Network), beta), RnaExpression(TumorNormalExp(TumorExpression, NormalExpression), Network)))
    }
    else if(identical(TumorExpression, FALSE)){
      print("Expression data are not provided, Only using Mutation data")
      RankResult = RankSum(Relevance(diffusedMatrix(normalizedMatrix(Network), beta), Mutation(SNP, Network)))
    }
    else{
      print("Using Mutation data and Expression data.")
      RankResult = RankSum(MinRelevance(Relevance(diffusedMatrix(normalizedMatrix(Network), beta), Mutation(SNP, Network)), Relevance(diffusedMatrix(normalizedMatrix(Network), beta), RnaExpression(TumorNormalExp(TumorExpression, NormalExpression), Network))))  
    }
  }
  else{
    print("Using your own network.")
    print("Canculate the restart probability of beta and this step needs more time.")
    beta = chooseBeta()
    print(paste0("the beta is ", beta))
    if (identical(SNP, FALSE) & identical(TumorExpression, FALSE)){
      print("You should provide Mutation data or Expression data.")
    }
    else if(identical(SNP, FALSE)){
      print("Mutation data are not provided, Only using Expression data")
      RankResult = RankSum(Relevance(diffusedMatrix(normalizedMatrix(Network), beta), RnaExpression(TumorNormalExp(TumorExpression, NormalExpression), Network)))
    }
    else if(identical(TumorExpression, FALSE)){
      print("Expression data are not provided, Only using Mutation data")
      RankResult = RankSum(Relevance(diffusedMatrix(normalizedMatrix(Network), beta), Mutation(SNP, Network)))
    }
    else{
      print("Using Mutation data and Expression data.")
      RankResult = RankSum(MinRelevance(Relevance(diffusedMatrix(normalizedMatrix(Network), beta), Mutation(SNP, Network)), Relevance(diffusedMatrix(normalizedMatrix(Network), beta), RnaExpression(TumorNormalExp(TumorExpression, NormalExpression), Network))))  
    }
  }
  
  if(identical(SNP, FALSE)){
    MinNetRankResult = matrix(0, length(RankResult), 3)
    rownames(MinNetRankResult) = names(RankResult)
    colnames(MinNetRankResult) = c("Gene", "Score", "KnownCancerGenes")
    MinNetRankResult[, 1] = names(RankResult)
    MinNetRankResult[, 2] = RankResult
    #knonwn cancer genes
    MinNetRankResult[intersect(rownames(MinNetRankResult), CGC[,1]), 3] = "CGC"
    MinNetRankResult[setdiff(rownames(MinNetRankResult), CGC[,1]), 3] = "-"
  }
  else {
    MinNetRankResult = matrix(0, length(RankResult), 4)
    rownames(MinNetRankResult) = names(RankResult)
    colnames(MinNetRankResult) = c("Gene", "Score", "KnownCancerGenes", "Freq")
    MinNetRankResult[, 1] = names(RankResult)
    MinNetRankResult[, 2] = RankResult
    #knonwn cancer genes
    MinNetRankResult[intersect(rownames(MinNetRankResult), CGC[,1]), 3] = "CGC"
    MinNetRankResult[setdiff(rownames(MinNetRankResult), CGC[,1]), 3] = "-"
    # mutation frequency
    MinNetRankResult[intersect(rownames(MinNetRankResult), rownames(SNP)), 4] =  apply(SNP,1,sum)[intersect(rownames(MinNetRankResult), rownames(SNP))]
    MinNetRankResult[,4] = as.numeric(MinNetRankResult[,4])/dim(SNP)[2]
  }
  return(MinNetRankResult)
}