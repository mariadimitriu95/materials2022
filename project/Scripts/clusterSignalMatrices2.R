#' clusterSignalMatrices
#'
#' @param ml A list of signal matrices, as produced by \code{\link{signal2Matrix}}.
#' @param k The number of clusters to generate
#' @param scaleRows Logical; whether to scale rows for clustering
#' @param scaleCols Logical; whether to scale columns (i.e. signals/samples)
#' @param assay Assay to use (ignored unless `ml` is an ESE object)
#' @param use What values to use for clustering. By default, uses 
#'   \code{\link[EnrichedHeatmap]{enriched_score}}. Other options are 'full' 
#'   (uses the full signal for clustering), 'max' (uses the maximum value in 
#'   the region), or 'center' (use the value at the center of the region).
#' @param by Optional factor/character/integer vector of the same length as 
#'   `ml`. When scaling rows, this can be used to indicate which rows should be
#'   scaled together.
#' @param trim Values to trim (applied individidually for each signal matrix)
#' @param nstart Number of starts for k-means clustering
#' @param ... Passed to `kmeans`
#'
#' @return If `k` is of length 1, a vector of cluster labels, corresponding to 
#'   the rows of `ml`. If `length(k)>1`, a list of two data.frames containing:
#'   1) the cluster labels at the different resolutions, and 2) the variance 
#'   explained by clusters at each resolution.
#' @export
#' @importFrom matrixStats rowMaxs rowVars
#' @importFrom stats kmeans
clusterSignalMatrices2 <- function(ml, k, scaleRows=FALSE, scaleCols=FALSE,
                                  use=c("enrich","full","max","center"),
                                  by=rep(1L,seq_along(ml)),
                                  assay=1L, trim=c(0.05,0.95), nstart=3, ...){
  if(is(ml, "SummarizedExperiment")) ml <- .ese2ml(ml, assay=assay)
  k <- unique(as.integer(k))
  stopifnot(all(k>1 & k<nrow(ml[[1]])))
  use <- match.arg(use)
  ml <- lapply(ml, FUN=function(x){
    q <- quantile(x, trim)
    x[x>q[2]] <- q[2]
    x[x<q[1]] <- q[1]
    x
  })
  if(scaleCols) ml <- lapply(ml, FUN=function(x) (x-mean(x))/sd(x))
  ml <- switch(use,
               full=ml,
               max=lapply(ml, FUN=matrixStats::rowMaxs),
               center=lapply(ml, FUN=function(x){
                 a <- attributes(x)
                 if(length(ti <- a$target_index)==0)
                   ti <- c(max(a$upstream_index),min(a$downstream_index))
                 rowMeans(x[,ti,drop=FALSE])
               }),
               enrich=lapply(ml, enriched_score))
  if(scaleRows){
    stopifnot(length(ml)==length(by))
    ml <- lapply(split(ml, by), FUN=function(x){
      m <- do.call(cbind, x)
      m <- m - rowMeans(m)
      sd <- sqrt(matrixStats::rowVars(m, na.rm=TRUE))
      sd[which(sd==0)] <- 1
      m <- m/sd
    })
  }
  m <- do.call(cbind, ml)
  res <- lapply(setNames(k,k), FUN=function(x){
    cl <- kmeans(dist(m), centers=k)
    ve <- round(100*sum(cl$betweenss)/sum(c(cl$withinss,cl$betweenss)))
    list(cl=factor(as.character(cl$cluster),as.character(seq_len(k))), ve=ve)
  })
  if(length(res)==1){
    message("  ~", res[[1]]$ve, "% of the variance explained by clusters")
    return(res[[1]]$cl)
  }
  cl <- do.call(cbind, lapply(res, FUN=function(x) as.data.frame(x$cl)))
  colnames(cl) <- k
  ve <- data.frame(k=k, varExplained=unlist(lapply(res, FUN=function(x) x$ve)))
  list(clusters=cl, varExplained=ve)
}
