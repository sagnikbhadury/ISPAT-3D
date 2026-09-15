# GP component validation: Vecchia approximation, scale sensitivity, and axial ablation.
suppressPackageStartupMessages({library(MASS);library(gpboost)})
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop('Usage: focused_gp_vecchia_validation.R output_dir')
outroot<-args[1]
set.seed(20260807)
N <- 300L; G <- 12L; n_rep <- 20L
strengths <- c(0,.5,1,1.5); kernels <- c('Matern','RBF')
lS <- .25; lZ <- .12
Omega <- diag(G)
for(g in seq_len(G)){h <- if(g==G) 1L else g+1L;Omega[g,h]<-Omega[h,g]<--.25}
Sigma <- solve(Omega); signal_var <- mean(diag(Sigma))
rho_true <- -Omega/sqrt(outer(diag(Omega),diag(Omega)));diag(rho_true)<-1
angles <- seq(0,2*pi,length.out=G+1L)[-1L]
L <- cbind(cos(angles),sin(angles))
kernel_matrix <- function(S,ls,lz,kernel) {
 dx <- outer(S[,1],S[,1],'-')/ls
 dy <- outer(S[,2],S[,2],'-')/ls
 dz <- outer(S[,3],S[,3],'-')/lz
 if(kernel=='Matern') {r<-sqrt(dx^2+dy^2+dz^2);K<-(1+sqrt(3)*r)*exp(-sqrt(3)*r)}
 else K<-exp(-.5*(dx^2+dy^2+dz^2))
 diag(K)<-1;K
}
knn_pairs <- function(S,k=8L){
 D<-as.matrix(dist(S));diag(D)<-Inf
 nbr<-t(apply(D,1,order))[,seq_len(k),drop=FALSE]
 list(i=rep(seq_len(nrow(S)),each=k),j=as.vector(t(nbr)),k=k)
}
moran <- function(R,pairs){
 R<-scale(R,center=TRUE,scale=FALSE)
 den<-colSums(R^2);val<-colSums(R[pairs$i,,drop=FALSE]*R[pairs$j,,drop=FALSE])/(pairs$k*den)
 mean(abs(val+1/(nrow(R)-1)),na.rm=TRUE)
}
pcor <- function(R){
 C<-cov(R);P<-solve(C+1e-8*mean(diag(C))*diag(ncol(C)))
 ans<--P/sqrt(outer(diag(P),diag(P)));diag(ans)<-1;ans
}
fixed_vecchia <- function(Y,coords,kernel,sd,ls,lz=NULL){
 if(sd==0) return(scale(Y,center=TRUE,scale=FALSE))
 d<-ncol(coords);fun<-if(kernel=='Matern') 'matern_ard' else 'gaussian_ard'
 pars<-c(signal_var,sd^2,rep(ls,2L),if(d==3L) lz)
 R<-matrix(0,nrow(Y),ncol(Y))
 for(g in seq_len(ncol(Y))){
  mod<-gpboost::fitGPModel(gp_coords=coords,cov_function=fun,
     cov_fct_shape=1.5,gp_approx='vecchia',num_neighbors=15L,
     num_parallel_threads=2L,y=Y[,g],X=matrix(1,nrow(Y),1),
     params=list(init_cov_pars=pars,estimate_cov_par_index=rep(0L,length(pars)),maxit=1L))
  pred<-predict(mod,gp_coords_pred=coords,X_pred=matrix(1,nrow(Y),1),
                predict_response=FALSE)[['mu']]
  R[,g]<-Y[,g]-pred
 }
 R
}
rows<-list();i<-1L
for(kernel in kernels) for(sd in strengths) for(rep in seq_len(n_rep)){
 S<-cbind(runif(N),runif(N),runif(N,0,.4))
 K<-kernel_matrix(S,lS,lZ,kernel)
 H<-t(chol(K+1e-8*diag(N)))%*%matrix(rnorm(N*2L),N,2L)
 B<-sd*H%*%t(L)
 Z<-MASS::mvrnorm(N,mu=rep(0,G),Sigma=Sigma)
 Y<-Z+B
 variants<-list(
   'Matched 3D Vecchia'=fixed_vecchia(Y,S,kernel,sd,lS,lZ),
   'Short scale 3D Vecchia'=fixed_vecchia(Y,S,kernel,sd,lS/2,lZ/2),
   'Planar 2D Vecchia'=fixed_vecchia(Y,S[,1:2,drop=FALSE],kernel,sd,lS),
   'No GP'=scale(Y,center=TRUE,scale=FALSE))
 pairs<-knn_pairs(S); upper<-upper.tri(rho_true);null<-upper & abs(rho_true)<1e-12
 Zc<-scale(Z,center=TRUE,scale=FALSE)
 for(method in names(variants)){
  R<-variants[[method]];rhat<-pcor(R)
  rows[[i]]<-data.frame(Kernel=kernel,Spatial_SD=sd,Replicate=rep,
    Method=method,Mean_abs_Moran=moran(R,pairs),
    Partial_correlation_RMSE=sqrt(mean((rhat[upper]-rho_true[upper])^2)),
    Null_edge_FPR=mean(abs(rhat[null])>.10),
    Signal_residual_RMSE=sqrt(mean((R-Zc)^2)))
  i<-i+1L
 }
 if(rep%%5L==0L) cat('DONE',kernel,sd,'rep',rep,'\n')
}
out<-data.table::rbindlist(rows)
dir.create(outroot,recursive=TRUE,showWarnings=FALSE)
data.table::fwrite(out,file.path(outroot,'focused_gp_vecchia_validation.csv'))
cat('COMPLETE FOCUSED GP VALIDATION\n')
