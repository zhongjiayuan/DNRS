#' MarkRank
markrank_revise <- function(dataset, matrix_adj, label, normal_num, alpha=0.85, eps=1e-8, E_value=NULL)
{
  
  
  m <- nrow(dataset)											# Sample number.
  n <- ncol(dataset)											# Gene number.
  if (length(label) != m) {
    stop("The size of label must accord with the sample number in dataset.")
  }
  
  M<-ncol(dataset)
  A<-matrix(0, nrow = M, ncol = M)
  net_dim<-dim(matrix_adj)
  for (i in 1:net_dim[1]){
    for (j in 2:net_dim[2]){
      value<-matrix_adj[i,j]
      center<-matrix_adj[i,1]
      nei<-matrix_adj[i,j]
      if (is.na(as.vector(value))){
        A[center,nei]<-0
        }else{
          label <- as.matrix(as.numeric(label))
          MI_center <- mpmi::mminjk(dataset[,center], label, level=0L, na.rm=FALSE)
          dataset_tmp_nei <- (dataset[,center] + dataset[,nei])/sqrt(2)
          MI_add_nei <- mpmi::mminjk(dataset_tmp_nei, label, level=0L, na.rm=FALSE)
          diff_MI_center_nei<-MI_add_nei-MI_center
          if (diff_MI_center_nei>0){
            A[center,nei]<-diff_MI_center_nei
          }else{
            A[center,nei]<-0
          }
          
        }
      }
  }
  
  diag(A) <- 0
  d<-rowSums(A)
  d[d==0]=1
  d[d!=1]=0
  degs <- rowSums(A)
  A <- t(A)
  for (i in 1:ncol(A)){
  if (i %% 1000 == 0) print(i)
  if (degs[i] != 0) A[,i] <- A[,i]/degs[i]
  }
  
  AA <-A
  
  if (class(E_value) == "NULL"){
    SD<-NULL
      for (i in 1:n){
        dat=dataset[,i]
        SD[i]<-stats::sd(dat[(normal_num+1):length(label)])
      }
    E_value <- abs(SD)/sum(SD)	
  }
  
  library(Hmisc)
  R1 <- 1
  R2 <- E_value
  tm <- 1  # Iteration steps.
  while (sum(abs(R1-R2)) >= eps){
    tm <- tm + 1
    R1 <- R2
    R2 <- alpha*AA%*%R1 +as.vector((alpha*(d%*%R1))%*%E_value)+(1-alpha)*E_value
    R2<-impute(R2, mean)
    
  }
  R <- as.vector(R2)/(mean(R2))
  names(R) <- colnames(dataset)
  
  return(list(score = R,
              steps = tm,
              NET2  = AA,
              initial_pars = list(alpha = alpha, eps = eps)
  ))
}

