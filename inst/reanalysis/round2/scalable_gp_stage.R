# Scalable anisotropic GP stage shared by the full 3D reanalysis and planar baseline.
suppressPackageStartupMessages({
  library(data.table)
  library(gpboost)
  library(jsonlite)
})
ZONE_NAMES <- c("Very Low","Low","Intermediate","High","Very High")
GP_NEIGHBORS <- 15L
GP_ANCHOR_MAX <- 5000L
GP_MAXIT <- 30L
GP_THREADS <- 4L

allocation <- function(counts, budget) {
  if (budget >= sum(counts)) return(as.integer(counts))
  raw <- as.double(budget) * as.double(counts) / sum(counts)
  n <- floor(raw)
  missing <- budget - sum(n)
  if (missing > 0L) {
    add <- order(raw-n,decreasing=TRUE)[seq_len(missing)]
    n[add] <- n[add]+1L
  }
  as.integer(n)
}
spatial_pick <- function(index, x, y, section, budget) {
  if (length(index) <= budget) return(index)
  groups <- split(index,section[index])
  nsec <- allocation(lengths(groups),budget)
  out <- integer()
  for (k in seq_along(groups)) {
    members <- groups[[k]]
    take <- nsec[k]
    if (take == 0L) next
    pool <- if (length(members) > 3L*take) sample(members,3L*take) else members
    xx <- x[pool]; yy <- y[pool]
    bx <- pmin(10L,1L+floor(10*(xx-min(xx))/(max(xx)-min(xx)+1e-9)))
    by <- pmin(10L,1L+floor(10*(yy-min(yy))/(max(yy)-min(yy)+1e-9)))
    bins <- split(pool,paste(bx,by,sep="_"))
    first <- vapply(bins,function(v) sample(v,1L),integer(1))
    if (length(first)>take) first <- sample(first,take)
    if (length(first)<take) {
      left <- setdiff(pool,first)
      first <- c(first,sample(left,min(take-length(first),length(left))))
    }
    out <- c(out,first)
  }
  if (length(out)>budget) out<-sample(out,budget)
  if (length(out)<budget) out<-c(out,sample(setdiff(index,out),budget-length(out)))
  out
}
gp_adjust_one <- function(y, coords, anchor, seed) {
  if (sd(y)<1e-9) {
    return(list(residual=y-mean(y),pars=rep(NA_real_,5),beta=mean(y),
                status="constant",fit_seconds=0,pred_seconds=0))
  }
  s_anchor <- coords[anchor,,drop=FALSE]
  y_anchor <- y[anchor]
  if (any(!is.finite(y_anchor))) stop("Non-finite anchor response")
  if (sd(y_anchor)<1e-9) {
    return(list(residual=y-mean(y),pars=rep(NA_real_,5),beta=mean(y),
                status="constant_anchor",fit_seconds=0,pred_seconds=0))
  }
  status <- "ok"
  fit_model <- function(params) gpboost::fitGPModel(
    gp_coords=s_anchor,cov_function="matern_ard",cov_fct_shape=1.5,
    gp_approx="vecchia",num_neighbors=GP_NEIGHBORS,
    num_parallel_threads=GP_THREADS,seed=as.integer(seed),
    y=y_anchor,X=matrix(1,length(y_anchor),1),params=params)
  tfit <- system.time({
    model <- tryCatch(
      fit_model(list(maxit=GP_MAXIT,delta_rel_conv=1e-3)),
      error=function(e) {
        status <<- "retry_init"
        v <- max(var(y_anchor),1e-4)
        ranges <- apply(s_anchor,2,function(a) max(diff(range(a))/3,1e-3))
        fit_model(list(maxit=GP_MAXIT,delta_rel_conv=1e-3,
                       init_cov_pars=c(max(0.2*v,1e-4),max(0.8*v,1e-4),ranges)))
      })
  })
  tpred <- system.time({
    fitted <- predict(
      model,gp_coords_pred=coords,X_pred=matrix(1,nrow(coords),1),
      predict_response=FALSE)[["mu"]]
  })
  if (length(fitted)!=length(y) || any(!is.finite(fitted)))
    stop("Invalid GP prediction length or values")
  list(residual=as.numeric(y-fitted),pars=as.numeric(model$get_cov_pars()),
       beta=as.numeric(model$get_coef()[1]),status=status,
       fit_seconds=tfit[["elapsed"]],pred_seconds=tpred[["elapsed"]])
}
make_manifest <- function(labels, zones, cov_dir, counts, filename) {
  records <- lapply(seq_along(zones),function(q) {
    tag <- gsub(" ","_",zones[q])
    list(name=zones[q],n=as.integer(counts[q]),
         cov=normalizePath(file.path(cov_dir,paste0("cov_",tag,".csv")),
                           winslash="/",mustWork=TRUE))
  })
  jsonlite::write_json(list(cell_types=labels,zones=records),
                       filename,auto_unbox=TRUE,pretty=TRUE)
}
