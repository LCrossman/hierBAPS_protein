preproc_alignment <-
function(snp.matrix){
  if(class(snp.matrix)!="matrix") stop("snp.matrix is not a valid matrix")

  n.seq <- nrow(snp.matrix)
prior <- apply(snp.matrix, 2, function(x) table(factor(x, levels = c("A", "R", "N", "D", "B", "C", "E", "Q", "Z", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")),exclude=c("-",NA)))
keep <- colSums(prior>0)>1
  
  if(sum(keep)==0){
    #all columns are conserved
    return(NA)
  }
  
  snp.matrix <- snp.matrix[, keep, drop=FALSE]
  prior <- prior[, keep, drop=FALSE]

  #finally generate a matrix for the prior nt values
  prior <- 1*(prior>0)
  prior <- t(t(prior)/colSums(prior))
orig.dist <- as.matrix(ape::dist.aa(ape::as.AAbin(snp.matrix), pairwise.deletion=TRUE))
return(list(
    n.seq = n.seq,
    dist = orig.dist,
    seq.inds = 1:n.seq,
    prior = prior,
    data = snp.matrix
  ))
}
