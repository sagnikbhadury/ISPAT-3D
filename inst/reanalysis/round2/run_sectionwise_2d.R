# Matched section-wise 2D ISPAT architecture with scalable planar Mat?rn GP.
# Uses precisely the cells selected for the middle 3D tier and the same
# log KDE transform, Vecchia approximation and downstream MSFA likelihood.
source('scalable_gp_stage.R')
options(warn=1)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=5L) stop('Usage: run_sectionwise_2d.R crc|bc budget input_kde_csv selected_ids_rds output_dir')
dataset <- args[1]; budget <- as.integer(args[2])
specs <- list(
 crc=list(data='INPUT_KDE_CSV',
  zcol='Z_um',section_col='section_number',labels=c('Tumor','CD8_T','CD4_T','Treg','T_cell','B_cell','Macrophage','Other_immune','Stroma')),
 bc=list(data='INPUT_KDE_CSV',
  zcol='Z',section_col=NULL,labels=c('Tumor_HER2pos','Tumor_basal','Tumor_luminal','Tumor_other','Endothelial','Macrophage','CAF','Myoepithelial','CD8_T_cell','CD4_T_cell','B_cell','Plasma_cell')))
if(!dataset %in% names(specs)) stop('Unknown dataset')
spec <- specs[[dataset]]
spec$data <- args[3]
selected <- readRDS(args[4])
out <- args[5]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
needed <- unique(c('X','Y',spec$zcol,spec$section_col,paste0('kde_',spec$labels)))
needed <- needed[!is.na(needed)]
cat('Reading source',dataset,'\n')
dat <- data.table::fread(spec$data,select=needed)
section <- if(is.null(spec$section_col)) as.integer(floor(dat[[spec$zcol]]/2)) else as.integer(dat[[spec$section_col]])
set.seed(2026)
gp2d <- function(y,coords,anchor,seed) {
 n <- length(y)
 if(n<10L || sd(y)<1e-9) return(list(residual=y-mean(y),status=if(n<10L) 'small_section' else 'constant',fit_s=0,pred_s=0))
 ya <- y[anchor]
 if(sd(ya)<1e-9) return(list(residual=y-mean(y),status='constant_anchor',fit_s=0,pred_s=0))
 status <- 'ok'
 runfit <- function(params) gpboost::fitGPModel(
  gp_coords=coords[anchor,,drop=FALSE],cov_function='matern_ard',cov_fct_shape=1.5,
  gp_approx='vecchia',num_neighbors=min(15L,length(anchor)-1L),
  
  num_parallel_threads=2L,seed=as.integer(seed),y=ya,
  X=matrix(1,length(ya),1),params=params)
 tf <- system.time({
  model <- tryCatch(runfit(list(maxit=30L,delta_rel_conv=1e-3)),
   error=function(e) {
    status <<- 'retry_init'
    v <- max(var(ya),1e-4)
    ran <- apply(coords[anchor,,drop=FALSE],2,function(x) max(diff(range(x))/3,1e-3))
    runfit(list(maxit=30L,delta_rel_conv=1e-3,init_cov_pars=c(.2*v,.8*v,ran)))
   })
 })
 tp <- system.time(pred <- predict(model,gp_coords_pred=coords,
                           X_pred=matrix(1,n,1),predict_response=FALSE)[['mu']])
 if(any(!is.finite(pred))) stop('Invalid planar GP prediction')
 list(residual=as.numeric(y-pred),status=status,fit_s=tf[['elapsed']],pred_s=tp[['elapsed']])
}
residuals <- vector('list',5L);names(residuals)<-ZONE_NAMES
meta <- list();counts <- integer(5L)
for(q in seq_along(ZONE_NAMES)) {
 tag <- gsub(' ','_',ZONE_NAMES[q]);pick <- as.integer(selected[[ZONE_NAMES[q]]])
 groups <- split(seq_along(pick),section[pick]);groups <- groups[lengths(groups)>0L]
 Z <- matrix(0,length(pick),length(spec$labels),dimnames=list(NULL,spec$labels))
 counts[q]<-length(pick)
 cat('ZONE',ZONE_NAMES[q],'n',length(pick),'sections',length(groups),'\n')
 for(i in seq_along(groups)) {
  loc <- groups[[i]];rows <- pick[loc]; xy <- as.matrix(dat[rows,c('X','Y'),with=FALSE])
  anchors <- if(length(loc)<=300L) seq_along(loc) else spatial_pick(seq_along(loc),xy[,1],xy[,2],
                               rep(1L,length(loc)),min(5000L,max(300L,ceiling(.1*length(loc)))))
  for(g in seq_along(spec$labels)) {
   label <- spec$labels[g]
   y <- log1p(dat[[paste0('kde_',label)]][rows]*1e9)
   fit <- tryCatch(gp2d(y,xy,anchors,seed=2026L+1000L*q+100L*i+g),
    error=function(e) list(residual=y-mean(y),status=paste0('fallback:',conditionMessage(e)),fit_s=NA_real_,pred_s=NA_real_))
   Z[loc,g]<-fit$residual
   meta[[length(meta)+1L]]<-data.frame(dataset=dataset,budget=budget,
    zone=ZONE_NAMES[q],section=names(groups)[i],cell_type=label,n=length(loc),
    anchors=length(anchors),status=fit$status,fit_s=fit$fit_s,pred_s=fit$pred_s)
  }
  if(i%%25L==0L || i==length(groups)) cat('  SECTIONS',i,'/',length(groups),'\n')
 }
 residuals[[q]]<-Z
 write.csv(cov(Z),file.path(out,paste0('cov_',tag,'.csv')))
}
saveRDS(residuals,file.path(out,'gp_residuals.rds'),compress=TRUE)
data.table::fwrite(data.table::rbindlist(meta),file.path(out,'gp_hyperparameters.csv'))
make_manifest(spec$labels,ZONE_NAMES,out,counts,file.path(out,'manifest.json'))
statuses <- table(data.table::rbindlist(meta)[['status']])
jsonlite::write_json(list(dataset=dataset,budget=budget,source_3d_selection=paste0('tier_',budget),
                         method='Matched section-wise planar GP Mat?rn ARD with Vecchia anchor likelihood and prediction',
                         sampled_counts=setNames(as.list(as.integer(counts)),ZONE_NAMES),
                         section_fits=length(meta),statuses=as.list(statuses)),
                    file.path(out,'provenance.json'),auto_unbox=TRUE,pretty=TRUE)
cat('COMPLETE 2D',dataset,budget,'\n')
