#Title: KNGP Algorithm
#Author: Chad Kimmel

#Summary: The following is the code for the KNGP Algorithm.  The user must input 4 sources of input.  
#Go to bottom of this script for an explanation of the inputs needed by the user and the main function call.  

require(foreach)
require(pROC)

#The main KGNPAlgorithm Call
KNGPAlgorithm <- function(LinkKnowledgeMatrix,NodeKnowledgeVector,NodeNamesVector,rootNodes)
{
  #Obtain the transition probability matrix. 
  TransitionProbMatrix <- getTransitionProbMatrix(LinkKnowledgeMatrix,apply(LinkKnowledgeMatrix,1,sum))
  rm(LinkKnowledgeMatrix); gc()
  
  #Obtain the root node vector index
  rootNodeVectorIndex <- foreach(i=1:length(rootNodes),.combine=c) %do% {which(NodeNamesVector==rootNodes[i])}
  
  #List of f values to test.  
  f_vector <- c(0,0.5,1,2,5,10,50,100,1000,10000)
  
  #Find the best_f
  best_f <- find_best_f(TransitionProbMatrix,NodeKnowledgeVector,rootNodeVectorIndex,f_vector)
  
  #Compute the prior prob vector using best_f
  priorProbVector <- getPriorProbVector(rootNodeVectorIndex,NodeKnowledgeVector,best_f)
  
  #Obtain the posterior prob vector.  
  posteriorProbVector <- inference(TransitionProbMatrix,matrix(priorProbVector,nrow=length(priorProbVector),ncol=1))
  
  #Pass back data frame
  result <- data.frame(cbind(name=NodeNamesVector,postProb=posteriorProbVector))
  names(result) <- c("name","postProb")
  result$name <- as.character(result$name)
  result$postProb <- as.numeric(as.character(result$postProb))
  result$name <- ifelse(result$name %in% rootNodes,paste0(result$name,"*"),result$name)
  return(list(KNGPResult=result,best_f=best_f))
  
}

#The find_best_f procedure.  
find_best_f <- function(TransitionProbMatrix,NodeKnowledgeVector,rootNodeVectorIndex,f_vector)
{
  best_f <- NULL
  best_AUC <- -100000
  
  nonRootNodes <- setdiff(1:dim(TransitionProbMatrix)[1],rootNodeVectorIndex) #Get the vector of all non-root nodes
  
  for (f in f_vector)
  {
    f_result <- foreach(i=1:length(rootNodeVectorIndex), .combine=rbind) %do%
    {
      rootIndex <- rootNodeVectorIndex[i] #The root index for the given loop
      proteinsToRank <- c(rootIndex, sample(nonRootNodes,getSampleSize(length(nonRootNodes))) ) #Create the list of proteins to rank.
      
      newRootNodes <- setdiff(rootNodeVectorIndex,rootIndex) #The new root nodes minus the rootIndex
      
      priorProbVector <- getPriorProbVector(newRootNodes,NodeKnowledgeVector,f) #Compute the prior prob vector
      posteriorProbVector <- inference(TransitionProbMatrix,matrix(priorProbVector,nrow=length(priorProbVector),ncol=1)) #Compute the posterior prob vector.
      data.frame(cbind(score=posteriorProbVector[,1][proteinsToRank],actual=c(1,rep(0,length(proteinsToRank)-1))))
    }
    auc <- as.numeric(auc(f_result$actual,f_result$score))
    if (auc>best_AUC) {best_AUC=auc;best_f=f}

  }
  
  return(best_f)
         
}

inference <- function(TransitionProbMatrix,priorProbMatrix)
{
  B = 0.5
  threshold = 0.00001
  delta = 1
  Po_prev = matrix(rep(0,dim(priorProbMatrix)[1]),nrow=dim(priorProbMatrix)[1],ncol=1)
  i = 1
  
  while (delta>threshold)
  {
    Po_curr <- ((1-B)*(TransitionProbMatrix %*% Po_prev)) + (B*priorProbMatrix)
    delta <- abs(sum(Po_curr) - sum(Po_prev))
    Po_prev <- Po_curr
#     print(paste("Iteration",i,"Delta",delta))
#     i <- i + 1
  }
  
  return(Po_prev)
}



#Get the transition probability matrix
getTransitionProbMatrix <- function(TransitionProbMatrix,sumMatrix)
{
  
  for (c in 1:dim(TransitionProbMatrix)[1])
  {
    if (sumMatrix[c]!=0) TransitionProbMatrix[,c] <- TransitionProbMatrix[,c]/sumMatrix[c] else
      TransitionProbMatrix[,c] <- rep(0,dim(TransitionProbMatrix)[1])  
  }
  
  return(TransitionProbMatrix)
}

getPriorProbVector <- function(rootNodes, NodeKnowledgeVector, f)
{
  newPriorValueVector <- foreach(s=1:length(NodeKnowledgeVector), .combine=c) %do% 
  {
    if (s %in% rootNodes) NodeKnowledgeVector[s]*f else
      NodeKnowledgeVector[s]
  }
  return(newPriorValueVector/sum(newPriorValueVector))
}

#Get the sample size to pull from.
getSampleSize <- function(size)
{
  if (size>=99) return(99) else
    return(size)
}

# #Define initial paramete values.  The 4 input parameters are the followin:
# # 1.  The n by n dimenstional link knowledge matrix
# # 2.  The n dimensional vector of node knowledge values
# # 3.  The n dimensional vector of node names
# # 4.  The vector of node names which are root nodes.  
# LinkKnowledgeMatrix <- matrix(data=c(1.0,0.3,0.3,0.1,0.3,0.3,1.0,0.6,0.2,0.3,0.3,0.6,1.0,0.1,0.4,0.1,0.2,0.1,1.0,0.2,0.3,0.3,0.4,0.2,1.0),nrow=5,ncol=5)
# NodeKnowledgeVector <- c(20,40,20,20,40)
# NodeNamesVector <- c("a","b","c","d","e")
# rootNodes <- c("a","b")
# 
# #Read in values from input file.  
# LinkKnowledgeCSV <- read.csv("C:\\PHDProject-GeneIdent\\Documents\\TechnicalPaper\\KNGPExampleInput\\LinkKnowledgeMatrix.csv",header=FALSE)
# LinkKnowledgeMatrix <- as.matrix(LinkKnowledgeCSV)
# 
# NodeKnowledgeCSV <- read.csv("C:\\PHDProject-GeneIdent\\Documents\\TechnicalPaper\\KNGPExampleInput\\NodeKnowledgeVector.csv",header=FALSE)
# NodeKnowledgeVector <- NodeKnowledgeCSV[,c(1)]
# 
# NodeNamesCSV <- read.csv("C:\\PHDProject-GeneIdent\\Documents\\TechnicalPaper\\KNGPExampleInput\\CandidateNodeSet.csv",header=FALSE)
# NodeNamesVector <- as.character(NodeNamesCSV[,c(1)])
# 
# RootNodesCSV <- read.csv("C:\\PHDProject-GeneIdent\\Documents\\TechnicalPaper\\KNGPExampleInput\\RootNodeSet.csv",header=FALSE)
# rootNodes <- as.character(RootNodesCSV[,c(1)])
# # 
##Dissertation Input
# NodeNamesCSV <- read.csv("C:\\PHDProject-GeneIdent\\Documents\\TechnicalPaper\\DissertationInput\\NodeNameVector.csv",header=FALSE)
# NodeNamesVector <- as.character(NodeNamesCSV[,c(1)])
# 
# NodeKnowledgeCSV <- read.csv("C:\\PHDProject-GeneIdent\\Documents\\TechnicalPaper\\DissertationInput\\NodeKnowledgeVector.csv",header=FALSE)
# NodeKnowledgeVector <- NodeKnowledgeCSV[,c(1)]
# 
# RootNodesCSV <- read.csv("C:\\PHDProject-GeneIdent\\Documents\\TechnicalPaper\\DissertationInput\\SeedListCSVFormat\\AlzheimerSeedList.csv",header=FALSE)
# rootNodes <- as.character(RootNodesCSV[,c(1)])
# 
# # # #Pass to KNGP procedure and produce output of posterior vaues ordered by Posterior Value.  
# KNGPResult <- KNGPAlgorithm(LinkKnowledgeMatrix,NodeKnowledgeVector,NodeNamesVector,rootNodes)$KNGPResult
# print(KNGPResult[order(KNGPResult$postProb,decreasing = TRUE),])

