reallocate_units_4 <- function(snp.object, partition, threshold=1e-5,
                               min.clust.size=20, split=FALSE,
                               n.cores=1){
  #some checks
  if (ncol(snp.object$prior)!=ncol(snp.object$data)) stop("ncol mismatch bwtn prior and data!")
  if (length(partition)!=nrow(snp.object$data)) stop("mismatch bwtn partition and data!")
  if (!(class(partition)=="integer")) stop("step 4 -> partition is not an integer vector!")
  clusters <- unique(partition)
  big_clusters <- clusters[purrr::map_int(clusters, ~ sum(partition==.x)>min.clust.size)==1]
  is.improved <- FALSE
  max_ml <- calc_log_ml(snp.object, partition)
  for(c in big_clusters){
    index <- c(1:length(partition))[partition==c]
    d <- snp.object$dist[partition==c, partition==c]
    h <- stats::hclust(stats::as.dist(d), method = "complete")
    if(split){
      npops <- 2
    } else {
      npops <- min(20, floor(nrow(d)/5))
    }
    sub.partition <- stats::cutree(h, k=npops)
    sub.clusters <- unique(sub.partition)
    max_ml <- calc_log_ml(snp.object, partition)
    temp_mls <- parallel::mclapply(sub.clusters, function(sub.c){
      best.clust <- calc_change_in_ml(snp.object, partition, index[sub.partition==sub.c])
      temp.partition <- partition
      temp.partition[index[sub.partition==sub.c]] <- as.integer(best.clust)
      return(list(lml=calc_log_ml(snp.object, temp.partition),
                  best.clust=best.clust))
    }, mc.cores=n.cores)
    arg.max <- which.max(purrr::map_dbl(temp_mls, ~ .x$lml))
    if(temp_mls[[arg.max]]$lml > (max_ml+threshold)){
      partition[index[sub.partition==sub.clusters[[arg.max]]]] <- as.integer(temp_mls[[arg.max]]$best.clust)
      max_ml <- temp_mls[[arg.max]]$lml
      is.improved <- TRUE
    }
  }
  return(list(partition=partition, is.improved=is.improved, lml=max_ml))
}
