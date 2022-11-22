save_lml_logs <- function(hb.object, file){
  if (!is.character(file)) stop("Invalid file name. Must be a valid string.")
  if ((!is.list(hb.object) || !is.data.frame(hb.object$partition.df)
       ) || !is.list(hb.object$lml.list)) stop("Invalid hb.object!")
  sink(file)
  for(i in 1:length(hb.object$lml.list)){
    cat(names(hb.object$lml.list)[[i]], "\n")
    cat("cluster names:\t", paste(names(hb.object$lml.list[[i]]), collapse = "\t"), "\n")
    cat("cluster log marginal likelihood:\t", paste(hb.object$lml.list[[i]], collapse = "\t"), "\n")
  }
  sink()
}
