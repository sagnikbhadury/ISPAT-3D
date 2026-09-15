# Section-block sufficient statistics from fixed middle-tier GP residuals.
source('scalable_gp_stage.R')
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=6L) stop('Usage: make_section_stats.R crc|bc budget input_kde_csv selected_ids_rds gp_residuals_rds output_json')
dataset<-args[1];budget<-as.integer(args[2])
root<-dirname(args[4])
selected<-readRDS(args[4])
residuals<-readRDS(args[5])
sourcefile<-args[3]
column<-if(dataset=='crc') 'section_number' else 'Z'
dat<-data.table::fread(sourcefile,select=column)
section<-if(dataset=='crc') as.integer(dat[[column]]) else as.integer(floor(dat[[column]]/2))
rm(dat);gc(verbose=FALSE)
zones<-vector('list',5L)
for(q in seq_along(ZONE_NAMES)){
 name<-ZONE_NAMES[q]
 ids<-as.integer(selected[[name]]);Z<-residuals[[name]]
 if(nrow(Z)!=length(ids)) stop('Selection and residual row mismatch in ',name)
 ids_sec<-section[ids]
 groups<-split(seq_along(ids),ids_sec)
 stats<-lapply(names(groups),function(s) {
  loc<-groups[[s]];A<-Z[loc,,drop=FALSE]
  C<-crossprod(A)
  list(section=as.integer(s),n=nrow(A),sum=as.numeric(colSums(A)),
       cross=lapply(seq_len(ncol(A)),function(i) as.numeric(C[i,])))
 })
 zones[[q]]<-list(name=name,records=stats)
 cat('ZONE',name,'sections',length(stats),'n',sum(vapply(stats,function(x)x$n,integer(1))),'\n')
}
sections<-sort(unique(section[unlist(selected,use.names=FALSE)]))
jsonlite::write_json(list(dataset=dataset,budget=budget,sections=as.integer(sections),
                          labels=colnames(residuals[[1]]),zones=zones),
 args[6],
 pretty=FALSE,auto_unbox=TRUE,digits=NA)
cat('COMPLETE SECTION STATS',dataset,'\n')
