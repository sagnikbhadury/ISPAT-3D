suppressPackageStartupMessages({library(data.table);library(jsonlite)})
source('scalable_gp_stage.R')
args <- commandArgs(trailingOnly=TRUE)
G <- as.integer(args[1])
if(length(args)!=4L) stop('Usage: run_saved_simulations.R G replicates simulation_source_dir output_root')
reps <- as.integer(args[2])
configs <- rbind(c(300,300,300),c(250,250,250),c(250,300,300),
                 c(250,250,300),c(300,250,250),c(300,300,250))
srcdir <- args[3]
outroot <- file.path(args[4],paste0('G_',G))
dir.create(outroot,recursive=TRUE,showWarnings=FALSE)
for (kernel in c('Matern','RBF')) for (k in seq_len(nrow(configs))) {
  nc <- configs[k,]
  file <- file.path(srcdir,paste0('simData3D_N_1_',nc[1],'_N_2_',nc[2],
    '_N_3_',nc[3],'_cells_',G,'_repli_5_',kernel,'_3D_datagenerate.rds'))
  if (!file.exists(file)) stop('Missing simulation file ',file)
  generated <- readRDS(file)
  for (r in seq_len(min(reps,length(generated)))) {
    tag <- paste(kernel,paste(nc,collapse='_'),paste0('rep_',r),sep='__')
    od <- file.path(outroot,tag)
    dir.create(od,recursive=TRUE,showWarnings=FALSE)
    done <- file.path(od,'manifest.json')
    if (file.exists(done)) next
    sim <- generated[[r]]
    YY <- sim$Gene_expression_matrix
    SS <- sim$S_loc_mat
    labels <- paste0('Marker_',seq_len(G))
    zone_records <- vector('list',3L)
    statuses <- list()
    for (q in 1:3) {
      idx <- which(SS[,4]==q)
      coords <- as.matrix(SS[idx,1:3,drop=FALSE])
      Z <- matrix(0,length(idx),G,dimnames=list(NULL,labels))
      for (g in seq_len(G)) {
        gp <- tryCatch(gp_adjust_one(as.numeric(YY[g,idx]),coords,
                        seq_along(idx),seed=2026L+10000L*G+100L*k+10L*q+g),
                error=function(e) list(residual=as.numeric(YY[g,idx]-mean(YY[g,idx])),
                                       status=paste0('fallback:',conditionMessage(e))))
        Z[,g] <- gp$residual
        statuses[[length(statuses)+1L]] <- list(zone=q,marker=g,status=gp$status)
      }
      write.csv(cov(Z),file.path(od,paste0('cov_Zone_',q,'.csv')))
      truth <- sim$Sigma_c[[q]]
      write.csv(truth,file.path(od,paste0('truth_Zone_',q,'.csv')))
      zone_records[[q]] <- list(name=paste0('Zone_',q),n=length(idx),
                              cov=normalizePath(file.path(od,paste0('cov_Zone_',q,'.csv')),
                                                winslash='/',mustWork=TRUE))
    }
    write.csv(tcrossprod(sim$Phi),file.path(od,'truth_Shared.csv'))
    write_json(list(G=G,kernel=kernel,config=as.integer(nc),replicate=r,statuses=statuses),
               file.path(od,'provenance.json'),pretty=TRUE,auto_unbox=TRUE)
    write_json(list(cell_types=labels,zones=zone_records),done,pretty=TRUE,auto_unbox=TRUE)
    cat('DONE',G,kernel,paste(nc,collapse=','),'rep',r,'\n')
  }
}
