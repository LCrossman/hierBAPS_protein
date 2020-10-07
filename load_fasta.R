load_fasta <-
function(alignment_file){

   library(ape) 
   seqs <- read.FASTA(alignment_file, type="AA")
   seq_names <- labels(seqs)
   seqs <- as.character(as.matrix(seqs))
   rownames(seqs) <- seq_names

   seqs[is.na(seqs)] <- "-"
   conserved <- colSums(t(t(seqs)==seqs[1,]))==nrow(seqs)

   seqs <- seqs[, !conserved]

   is_singleton <- apply(seqs,2,function(x){
       tab <- table(x)
       return(x %in% names(tab)[tab==1])
       })

   seqs[is_singleton] <- "-"
   seqs[seqs=="X"] <- "-"

   snp_mat <- seqs
   return(snp_mat)

}