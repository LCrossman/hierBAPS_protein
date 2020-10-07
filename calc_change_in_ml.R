calc_change_in_ml <-
function(snp.object, partition, indexes){
  #some checks
  if (ncol(snp.object$prior)!=ncol(snp.object$data)) stop("ncol mismatch bwtn prior and data!")
  if (length(partition)!=nrow(snp.object$data)) stop("mismatch bwtn partition and data!")
  if (!(all(indexes %in% 1:length(partition)))) stop("indexes outside of partiton range!")
  if (!(class(partition)=="integer")) stop("partition is not an integer vector!")

  original_cluster <- unique(partition[indexes])
  if(length(original_cluster)!=1) stop("there was not a unique cluster in the index set!")

  #create temporary partition with indexes in a seperate cluster
  temp_partition <- partition
  temp_cluster <- max(partition)+1
  temp_partition[indexes] <- temp_cluster
      mA <- t(rowsum(1*(snp.object$data=="A"),temp_partition))
      mR <- t(rowsum(1*(snp.object$data=="R"),temp_partition))
      mN <- t(rowsum(1*(snp.object$data=="N"),temp_partition))
      mD <- t(rowsum(1*(snp.object$data=="D"),temp_partition))
      mB <- t(rowsum(1*(snp.object$data=="B"),temp_partition))
      mC <- t(rowsum(1*(snp.object$data=="C"),temp_partition))
      mE <- t(rowsum(1*(snp.object$data=="E"),temp_partition))
      mQ <- t(rowsum(1*(snp.object$data=="Q"),temp_partition))
      mZ <- t(rowsum(1*(snp.object$data=="Z"),temp_partition))
      mG <- t(rowsum(1*(snp.object$data=="G"),temp_partition))
      mH <- t(rowsum(1*(snp.object$data=="H"),temp_partition))
      mI <- t(rowsum(1*(snp.object$data=="I"),temp_partition))
      mL <- t(rowsum(1*(snp.object$data=="L"),temp_partition))
      mK <- t(rowsum(1*(snp.object$data=="K"),temp_partition))
      mM <- t(rowsum(1*(snp.object$data=="M"),temp_partition))
      mF <- t(rowsum(1*(snp.object$data=="F"),temp_partition))
      mP <- t(rowsum(1*(snp.object$data=="P"),temp_partition))
      mS <- t(rowsum(1*(snp.object$data=="S"),temp_partition))
      mT <- t(rowsum(1*(snp.object$data=="T"),temp_partition))
      mW <- t(rowsum(1*(snp.object$data=="W"),temp_partition))
      mY <- t(rowsum(1*(snp.object$data=="Y"),temp_partition))
      mV <- t(rowsum(1*(snp.object$data=="V"),temp_partition))

      mA <- mA[,colnames(mA)!=temp_cluster, drop=FALSE] + mA[, colnames(mA)==temp_cluster]
      mR <- mR[,colnames(mR)!=temp_cluster, drop=FALSE] + mR[, colnames(mR)==temp_cluster]
      mN <- mN[,colnames(mN)!=temp_cluster, drop=FALSE] + mN[, colnames(mN)==temp_cluster]
      mD <- mD[,colnames(mD)!=temp_cluster, drop=FALSE] + mD[, colnames(mD)==temp_cluster]
      mB <- mB[,colnames(mB)!=temp_cluster, drop=FALSE] + mB[, colnames(mB)==temp_cluster]
      mC <- mC[,colnames(mC)!=temp_cluster, drop=FALSE] + mC[, colnames(mC)==temp_cluster]
      mE <- mE[,colnames(mE)!=temp_cluster, drop=FALSE] + mE[, colnames(mE)==temp_cluster]
      mQ <- mQ[,colnames(mQ)!=temp_cluster, drop=FALSE] + mQ[, colnames(mQ)==temp_cluster]
      mZ <- mZ[,colnames(mZ)!=temp_cluster, drop=FALSE] + mZ[, colnames(mZ)==temp_cluster]
      mG <- mG[,colnames(mG)!=temp_cluster, drop=FALSE] + mG[, colnames(mG)==temp_cluster]
      mH <- mH[,colnames(mH)!=temp_cluster, drop=FALSE] + mH[, colnames(mH)==temp_cluster]
      mI <- mI[,colnames(mI)!=temp_cluster, drop=FALSE] + mI[, colnames(mI)==temp_cluster]
      mL <- mL[,colnames(mL)!=temp_cluster, drop=FALSE] + mL[, colnames(mL)==temp_cluster]
      mK <- mK[,colnames(mK)!=temp_cluster, drop=FALSE] + mK[, colnames(mK)==temp_cluster]
      mM <- mM[,colnames(mM)!=temp_cluster, drop=FALSE] + mM[, colnames(mM)==temp_cluster]
      mF <- mF[,colnames(mF)!=temp_cluster, drop=FALSE] + mF[, colnames(mF)==temp_cluster]
      mP <- mP[,colnames(mP)!=temp_cluster, drop=FALSE] + mP[, colnames(mP)==temp_cluster]
      mS <- mS[,colnames(mS)!=temp_cluster, drop=FALSE] + mS[, colnames(mS)==temp_cluster]
      mT <- mT[,colnames(mT)!=temp_cluster, drop=FALSE] + mT[, colnames(mT)==temp_cluster]
      mW <- mW[,colnames(mW)!=temp_cluster, drop=FALSE] + mW[, colnames(mW)==temp_cluster]
      mY <- mY[,colnames(mY)!=temp_cluster, drop=FALSE] + mY[, colnames(mY)==temp_cluster]
      mV <- mV[,colnames(mV)!=temp_cluster, drop=FALSE] + mV[, colnames(mV)==temp_cluster]

   prior <- snp.object$prior
   prior[prior==0] <- 1 #deal with zeros and resulting NAs
 
   #calculate log marginal likelihood
   term1 <- -lgamma(1 + mA+mR+mN+mD+mB+mC+mE+mQ+mZ+mG+mH+mI+mL+mK+mM+mF+mP+mS+mT+mW+mY+mV)
 term2 <- lgamma(prior["A", ] + mA) - lgamma(prior["A", ])
 term2 <- term2 + lgamma(prior["R", ] + mR) - lgamma(prior["R", ])
 term2 <- term2 + lgamma(prior["N", ] + mN) - lgamma(prior["N", ])
 term2 <- term2 + lgamma(prior["D", ] + mD) - lgamma(prior["D", ])
 term2 <- term2 + lgamma(prior["B", ] + mB) - lgamma(prior["B", ])
 term2 <- term2 + lgamma(prior["C", ] + mC) - lgamma(prior["C", ])
 term2 <- term2 + lgamma(prior["E", ] + mE) - lgamma(prior["E", ])
 term2 <- term2 + lgamma(prior["Q", ] + mQ) - lgamma(prior["Q", ])
 term2 <- term2 + lgamma(prior["Z", ] + mZ) - lgamma(prior["Z", ])
 term2 <- term2 + lgamma(prior["G", ] + mG) - lgamma(prior["G", ])
 term2 <- term2 + lgamma(prior["H", ] + mH) - lgamma(prior["H", ])
 term2 <- term2 + lgamma(prior["I", ] + mI) - lgamma(prior["I", ])
 term2 <- term2 + lgamma(prior["L", ] + mL) - lgamma(prior["L", ])
 term2 <- term2 + lgamma(prior["K", ] + mK) - lgamma(prior["K", ])
 term2 <- term2 + lgamma(prior["M", ] + mM) - lgamma(prior["M", ])
 term2 <- term2 + lgamma(prior["F", ] + mF) - lgamma(prior["F", ])
 term2 <- term2 + lgamma(prior["P", ] + mP) - lgamma(prior["P", ])
 term2 <- term2 + lgamma(prior["S", ] + mS) - lgamma(prior["S", ])
 term2 <- term2 + lgamma(prior["T", ] + mT) - lgamma(prior["T", ])
 term2 <- term2 + lgamma(prior["W", ] + mW) - lgamma(prior["W", ])
 term2 <- term2 + lgamma(prior["Y", ] + mY) - lgamma(prior["Y", ])
 term2 <- term2 + lgamma(prior["V", ] + mV) - lgamma(prior["V", ])

  new_columnwise_lml <- colSums(term1 + term2)
  names(new_columnwise_lml) <- colnames(mA)

  mA <- t(rowsum(1*(snp.object$data=="A"),partition))
      mR <- t(rowsum(1*(snp.object$data=="R"),partition))
      mN <- t(rowsum(1*(snp.object$data=="N"),partition))
      mD <- t(rowsum(1*(snp.object$data=="D"),partition))
      mB <- t(rowsum(1*(snp.object$data=="B"),partition))
      mC <- t(rowsum(1*(snp.object$data=="C"),partition))
      mE <- t(rowsum(1*(snp.object$data=="E"),partition))
      mQ <- t(rowsum(1*(snp.object$data=="Q"),partition))
      mZ <- t(rowsum(1*(snp.object$data=="Z"),partition))
      mG <- t(rowsum(1*(snp.object$data=="G"),partition))
      mH <- t(rowsum(1*(snp.object$data=="H"),partition))
      mI <- t(rowsum(1*(snp.object$data=="I"),partition))
      mL <- t(rowsum(1*(snp.object$data=="L"),partition))
      mK <- t(rowsum(1*(snp.object$data=="K"),partition))
      mM <- t(rowsum(1*(snp.object$data=="M"),partition))
      mF <- t(rowsum(1*(snp.object$data=="F"),partition))
      mP <- t(rowsum(1*(snp.object$data=="P"),partition))
      mS <- t(rowsum(1*(snp.object$data=="S"),partition))
      mT <- t(rowsum(1*(snp.object$data=="T"),partition))
      mW <- t(rowsum(1*(snp.object$data=="W"),partition))
      mY <- t(rowsum(1*(snp.object$data=="Y"),partition))
      mV <- t(rowsum(1*(snp.object$data=="V"),partition))


 term1 <- -lgamma(1 + mA+mR+mN+mD+mB+mC+mE+mQ+mZ+mG+mH+mI+mL+mK+mM+mF+mP+mS+mT+mW+mY+mV)
 term2 <- lgamma(prior["A", ] + mA) - lgamma(prior["A", ])
 term2 <- term2 + lgamma(prior["R", ] + mR) - lgamma(prior["R", ])
 term2 <- term2 + lgamma(prior["N", ] + mN) - lgamma(prior["N", ])
 term2 <- term2 + lgamma(prior["D", ] + mD) - lgamma(prior["D", ])
 term2 <- term2 + lgamma(prior["B", ] + mB) - lgamma(prior["B", ])
 term2 <- term2 + lgamma(prior["C", ] + mC) - lgamma(prior["C", ])
 term2 <- term2 + lgamma(prior["E", ] + mE) - lgamma(prior["E", ])
 term2 <- term2 + lgamma(prior["Q", ] + mQ) - lgamma(prior["Q", ])
 term2 <- term2 + lgamma(prior["Z", ] + mZ) - lgamma(prior["Z", ])
 term2 <- term2 + lgamma(prior["G", ] + mG) - lgamma(prior["G", ])
 term2 <- term2 + lgamma(prior["H", ] + mH) - lgamma(prior["H", ])
 term2 <- term2 + lgamma(prior["I", ] + mI) - lgamma(prior["I", ])
 term2 <- term2 + lgamma(prior["L", ] + mL) - lgamma(prior["L", ])
 term2 <- term2 + lgamma(prior["K", ] + mK) - lgamma(prior["K", ])
 term2 <- term2 + lgamma(prior["M", ] + mM) - lgamma(prior["M", ])
 term2 <- term2 + lgamma(prior["F", ] + mF) - lgamma(prior["F", ])
 term2 <- term2 + lgamma(prior["P", ] + mP) - lgamma(prior["P", ])
 term2 <- term2 + lgamma(prior["S", ] + mS) - lgamma(prior["S", ])
 term2 <- term2 + lgamma(prior["T", ] + mT) - lgamma(prior["T", ])
 term2 <- term2 + lgamma(prior["W", ] + mW) - lgamma(prior["W", ])
 term2 <- term2 + lgamma(prior["Y", ] + mY) - lgamma(prior["Y", ])
 term2 <- term2 + lgamma(prior["V", ] + mV) - lgamma(prior["V", ])

orignal_columnwise_lml <- colSums(term1 + term2)
  names(orignal_columnwise_lml) <- colnames(mA)
  new_columnwise_lml <- new_columnwise_lml[names(orignal_columnwise_lml)]

  diff_columnwise_lml <- orignal_columnwise_lml - new_columnwise_lml
  diff_columnwise_lml[original_cluster] <- Inf

  best_column <- names(orignal_columnwise_lml)[which.min(diff_columnwise_lml)]

  return(best_column)
}
