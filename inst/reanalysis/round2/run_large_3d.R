source("scalable_gp_stage.R")
options(warn=1)
ROOT <- getwd()
SPECS <- list(
  crc=list(
    data="INPUT_KDE_CSV",
    labels=c("Tumor","CD8_T","CD4_T","Treg","T_cell","B_cell","Macrophage","Other_immune","Stroma"),
    zcol="Z_um",section_col="section_number",zone_col="pathology_zone",
    budgets=c(50000L,100000L,150000L),budget_kind="per_zone"),
  bc=list(
    data="INPUT_KDE_CSV",
    labels=c("Tumor_HER2pos","Tumor_basal","Tumor_luminal","Tumor_other",
             "Endothelial","Macrophage","CAF","Myoepithelial",
             "CD8_T_cell","CD4_T_cell","B_cell","Plasma_cell"),
    zcol="Z",section_col=NULL,zone_col="pathology_zone",
    budgets=c(25000L,45000L,65000L),budget_kind="total")
)
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L || !args[1] %in% names(SPECS))
  stop("Usage: run_large_3d.R crc|bc input_kde_csv output_root")
dataset <- args[1]
spec <- SPECS[[dataset]]
spec$data <- args[2]
out_root <- file.path(args[3],dataset)
dir.create(out_root,recursive=TRUE,showWarnings=FALSE)
needed <- c("X","Y",spec$zcol,spec$zone_col,spec$section_col,
            paste0("kde_",spec$labels))
needed <- unique(needed[!is.na(needed)])
cat("Reading",dataset,"source",spec$data,"\n")
tm <- system.time(dat <- data.table::fread(spec$data,select=needed))
cat("Rows",nrow(dat),"read_seconds",tm[["elapsed"]],"\n")
zone <- dat[[spec$zone_col]]
if (is.character(zone)) zone <- match(zone,ZONE_NAMES)
zone <- as.integer(zone)
if (anyNA(zone) || any(!zone %in% 1:5)) stop("Invalid zone labels")
section <- if (is.null(spec$section_col))
  as.integer(floor(dat[[spec$zcol]]/2)) else
  as.integer(dat[[spec$section_col]])
if (anyNA(section)) stop("Invalid section IDs")
all_count <- tabulate(zone,5)
cat("Zone populations:",paste(ZONE_NAMES,all_count,sep="=",collapse=", "),"\n")
xy_x <- dat[["X"]]; xy_y <- dat[["Y"]]
for (tier in seq_along(spec$budgets)) {
  budget <- spec$budgets[tier]
  n_zone <- if (spec$budget_kind=="per_zone")
    rep(budget,5) else allocation(all_count,budget)
  tier_dir <- file.path(out_root,paste0("tier_",budget))
  dir.create(tier_dir,recursive=TRUE,showWarnings=FALSE)
  cat("TIER",budget,"counts",paste(n_zone,collapse=","),"\n")
  selected <- vector("list",5)
  residuals <- vector("list",5)
  names(selected) <- names(residuals) <- ZONE_NAMES
  meta <- list()
  set.seed(2026L+tier)
  for (q in 1:5) {
    tag <- gsub(" ","_",ZONE_NAMES[q])
    idx <- which(zone==q)
    if (n_zone[q]>length(idx)) stop("Budget exceeds available cells for ",ZONE_NAMES[q])
    pick <- spatial_pick(idx,xy_x,xy_y,section,n_zone[q])
    selected[[q]] <- pick
    coords <- as.matrix(dat[pick,c("X","Y",spec$zcol),with=FALSE])
    storage.mode(coords) <- "double"
    anchors <- spatial_pick(seq_along(pick),coords[,1],coords[,2],
                            section[pick],min(GP_ANCHOR_MAX,max(300L,as.integer(ceiling(0.10*length(pick))))))
    Z <- matrix(0,length(pick),length(spec$labels),
                dimnames=list(NULL,spec$labels))
    cat("  ZONE",ZONE_NAMES[q],"n",length(pick),"anchors",length(anchors),
        "sections",length(unique(section[pick])),"\n")
    for (g in seq_along(spec$labels)) {
      label <- spec$labels[g]
      y <- log1p(dat[[paste0("kde_",label)]][pick]*1e9)
      fit <- gp_adjust_one(y,coords,anchors,seed=2026L+100L*tier+10L*q+g)
      Z[,g] <- fit$residual
      row <- data.frame(dataset=dataset,tier=budget,zone=ZONE_NAMES[q],
                        cell_type=label,n=length(pick),anchors=length(anchors),
                        status=fit$status,fit_s=fit$fit_seconds,predict_s=fit$pred_seconds,
                        beta=fit$beta,residual_variance=var(fit$residual),
                        nugget=fit$pars[1],GP_variance=fit$pars[2],
                        range_X=fit$pars[3],range_Y=fit$pars[4],range_Z=fit$pars[5])
      meta[[length(meta)+1L]] <- row
      cat("    GP",label,"fit_s",round(fit$fit_seconds,2),
          "pred_s",round(fit$pred_seconds,2),
          "resvar",signif(var(fit$residual),3),"\n")
    }
    residuals[[q]] <- Z
    C <- cov(Z)
    write.csv(C,file.path(tier_dir,paste0("cov_",tag,".csv")))
    rm(idx,pick,coords,anchors,Z,C); gc(verbose=FALSE)
  }
  saveRDS(selected,file.path(tier_dir,"selected_ids.rds"),compress=TRUE)
  saveRDS(residuals,file.path(tier_dir,"gp_residuals.rds"),compress=TRUE)
  data.table::fwrite(data.table::rbindlist(meta),
                     file.path(tier_dir,"gp_hyperparameters.csv"))
  make_manifest(spec$labels,ZONE_NAMES,tier_dir,n_zone,
                file.path(tier_dir,"manifest.json"))
  jsonlite::write_json(list(dataset=dataset,tier=budget,budget_kind=spec$budget_kind,
                            sampled_counts=setNames(as.list(as.integer(n_zone)),ZONE_NAMES),
                            source_rows=nrow(dat),method="5,000 balanced anchor Mat?rn ARD GP with Vecchia likelihood"),
                       file.path(tier_dir,"run_provenance.json"),
                       auto_unbox=TRUE,pretty=TRUE)
  cat("COMPLETED TIER",budget,"\n")
}
cat("ALL",dataset,"GP TIERS COMPLETE\n")
