split_clusters_3 <- function(snp.object, partition, threshold=1e-5,
                             min.clust.size=20, n.cores=1){
  #At the moment this can't create new clusters. This is the same as in the original hierBAPS
  #but it might be worth allowing the creation of new clusters. TODO:Ask Jukka about it.
  return(reallocate_units_4(snp.object, partition, threshold=1e-5,
                            min.clust.size=min.clust.size, split=TRUE,
                            n.cores=n.cores))
}
