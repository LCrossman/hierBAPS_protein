move_units_1 <- function(snp.object, partition, threshold=1e-5,
                         frac.clust.searched=0.3,
                         min.clust.size=20,
                         n.cores=1){
  #some checks
  if (ncol(snp.object$prior)!=ncol(snp.object$data)) stop("ncol mismatch bwtn prior and data!")
  if (length(partition)!=nrow(snp.object$data)) stop("mismatch bwtn partition and data!")
  if (!(class(partition)=="integer")) stop("in move 1 -> partition is not an integer vector!")
  nsamples <- nrow(snp.object$data)
  clusters <- unique(partition)
  cluster_size <- purrr::map_int(clusters, ~ sum(partition==.x))
  names(cluster_size) <- clusters
  max_ml <- calc_log_ml(snp.object, partition)
  max_partition <- partition
  #identify the most divereged isolates of each cluster as candidates for moving.
  #as in hierBAPS we take the individuals that are most distant from each other
  indexes <- purrr::as_vector(purrr::imap(clusters[cluster_size>1], function(c, ind){
    if(cluster_size[ind] > min.clust.size){
      #big cluster so only move the most distant isolates
      index <- c(1:length(partition))[partition==c]
      d <- snp.object$dist[partition==c, partition==c]
      diag(d) <- NA
      mean.dist <- rowMeans(d, na.rm = TRUE)
      return(index[mean.dist > stats::quantile(mean.dist, probs=(1-frac.clust.searched))])
    } else {
      #small cluster so consider all isolates
      return(c(1:length(partition))[partition==c])
    }
  }))
  is.improved <- FALSE
  for (i in indexes){
    best.cluster <- calc_change_in_ml(snp.object, max_partition, i)
    temp_partition <- max_partition
    temp_partition[i] <- as.integer(best.cluster)
    temp_lml <- calc_log_ml(snp.object, temp_partition)
    #If improvment above threshold make the swap
    if(temp_lml > (max_ml+threshold)){
      max_ml <- temp_lml
      max_partition <- temp_partition
      is.improved <- TRUE
    }
  }
  return(list(partition=max_partition, is.improved=is.improved, lml=max_ml))
}
