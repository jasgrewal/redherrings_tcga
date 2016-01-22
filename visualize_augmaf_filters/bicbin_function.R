#Improve subsampling routine
score_C <- function(sub_d_m, sub_d_n,d_sum,p) { #d_orig is subsampled matrix, d_sum is sum from submatrix of collapsed sums, p is average proportion of ones
  dimprod=sub_d_m*sub_d_n*p #Expected ones in this submatrix?
  k=d_sum #How many ones actually in this submatrix?
  #Calculate score
  if((k >= dimprod) & (k <= (2*dimprod))) {
    return ((k-dimprod)^2)/(3*(sub_d_m+1)*(sub_d_n+1)*p)
  } else {
    if(k > (2*dimprod)){
      return ((k-dimprod)^2)/(sub_d_m*sub_d_n*(k + dimprod))
    } else{
      return (0) #Only defined for k >= mnp?
    }
  }
}

subsample_routine <- function(d, n=10,len=1,colrange=seq_len(ncol(d)), rowrange=seq_len(nrow(d))){#sampling n points with replacement (r flag), len=1 for col, 2 for row based
  set.seed(1)
  #Subsample both row and columns randomly (yields square matrix)
  #sub_d <- d[sample(seq_len(nrow(d)),size=n),sample(seq_len(ncol(d)),size=n)]
  if(len==1){ #If subselecting all rows for random columns
    sub_d <- d[,sample(seq_len(ncol(d)),size=n,replace=FALSE),drop=FALSE]
  } else{
    sub_d <- d[sample(seq_len(nrow(d)),size=n,replace=FALSE),,drop=FALSE]
  }
  return (remove_all_zeroed(sub_d))
}

update_mat <- function(update_df,ones_d,rc,score_add=0){
  if(rc == 1){
    #Subsampled by columns, so getting the arg max row values
    update_df=update_df[rownames(update_df) %in% rownames(ones_d),] + score_add
  } else{
    #Subsampled by rows, so getting the arg max col values
    update_df=update_df[,colnames(update_df) %in% rownames(ones_d)] + score_add
  }
  return(update_df)
}

get_ones<-function(mymat){
  return((length(mymat[mymat > 0]))/(length(mymat)))
}

remove_all_zeroed <- function(mymat){
  mymat=mymat[which(rowSums(mymat)>0),which(colSums(mymat)>0)]
  mymat1=mymat[,which(colSums(mymat) <= quantile(colSums(mymat),0.95))]
  mymat2=mymat1[which(rowSums(mymat1)>1),]
  return(mymat2)
}

normalize_mat <- function(mymat){
  temp=mymat;
  for(i in 1:ncol(mymat)){
    #Adjust each col's (all genes for patient i) score based on total non zero entries for patient i
    mymat[,i] = mymat[,i,drop=FALSE]/(colSums(mymat[,i,drop=FALSE])[[1]])
  }
  for(j in 1:nrow(mymat)){
    #Now multiply each variant entry based on the cohort wide occurrence counts
    mymat[j,] = mymat[j,,drop=FALSE] * dim(mymat[j,which(mymat[j,,drop=FALSE]>0),drop=FALSE])[2]
  }
  return(mymat)
}

trim_mat=normalize_mat(remove_all_zeroed(binary_matrix));
trim_mat_ordered=trim_mat[order(rowSums(trim_mat),decreasing=TRUE),]
trim_mat_ordered=trim_mat[,order(colSums(trim_mat),decreasing=TRUE)]
trim_mat_ordered_normalized=normalize(trim_mat_ordered)
library(scatterplot3d)
nc=ncol(trim_mat);nl=nrow(trim_mat)
scatterplot3d(rep(1:nl,nc),rep(1:nc,each=nl), as.vector(trim_mat_ordered),
              col.axis="blue",angle=40,
              col.grid="lightblue", main="", xlab="", ylab="", zlab="",
              pch=21, box=FALSE, cex.symbols=1,type="h",color="red",axis=FALSE)


#Iterative algorithm to find maximal C bicluster
search_max_C <- function(d, subsize=100) { 
  avg_ones= get_ones(d)
  my_blind_dat = d; my_blind_dat[] <- 0
  iter_num=0
  while((sum(my_blind_dat)==0 ) & (iter_num<11) ){
    iter_num=iter_num+1
    for(sub_by in 1:2){
      #maximize over a subset of cols, then a subset of rows
      #1 is cols, 2 is rows
      sub_d <- subsample_routine(d,n=subsize,len=sub_by) 
      #Get collapsed sums over subsampled rows/columns. Get 1 x 'x' vector
      if(sub_by==1){
        sub_d_sum = as.matrix(t(rowSums(sub_d))); sub_d=t(sub_d)
      } else{
        sub_d_sum = as.matrix(t(colSums(sub_d)))
      }
      #Reorder by decreasing values
      sub_d_sum<- as.matrix(t(sub_d_sum[,order(colnames(sub_d_sum),decreasing=TRUE)]));
      max_c_vector <- character(0)
      for (i in 1:length(sub_d_sum)) {
        c_score=score_C(sub_d_m=i, sub_d_n=dim(sub_d)[2],sum(sub_d_sum[1:i]),p=avg_ones)
        max_c_vector <- c(max_c_vector,c_score)
      }
      
      max_c_mat <- as.matrix(max_c_vector);rownames(max_c_mat)<-colnames(sub_d_sum)
      max_score = as.numeric(max(max_c_mat));
      max_c_mat=max_c_mat[which(max_c_mat==max(max_c_mat)),,drop=FALSE];
      my_blind_dat=update_mat(my_blind_dat,max_c_mat, rc=sub_by, score_add=max_score) #Update original matrix with bicluster highest
      if(sum(my_blind_dat)>0){print(dim(max_c_mat)); print(max_score); if(dim(max_c_mat)[1]<3){print(max_c_mat)}}
      
    }
    print(iter_num)
    if(sum(my_blind_dat)>0){ break}
    
  }
  return(sum(my_blind_dat))
}

for(subs in seq(50,913,by=50)){
  print("--------"); print(paste("SUBSET SIZE: ",subs)); print("--------")
  search_max_C(binary_matrix,subsize=subs)
}
