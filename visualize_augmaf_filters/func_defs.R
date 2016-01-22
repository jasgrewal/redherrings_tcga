makefreqdf <- function(mydf,distype) {
  t_mutc=unique(mydf[,c("Hugo_Symbol","mutation","Tumor_Sample_Barcode")]);
  tf_genc=count(t_mutc, c('Hugo_Symbol')); tf_mutc=count( t_mutc, c('Hugo_Symbol','Tumor_Sample_Barcode')); 
  tf_mutc$gene_count<-tf_genc[match(tf_mutc$Hugo_Symbol,tf_genc$Hugo_Symbol),2]; 
  tf_mutc$mut_fcount<-as.factor(tf_mutc$freq)
  tf_mutc<-tf_mutc[with(tf_mutc, order(-(gene_count),-(freq),Hugo_Symbol)),];tf_mutc$Hugo_Symbol<-factor(tf_mutc$Hugo_Symbol,levels=tf_mutc$Hugo_Symbol)
  tf_mutc$disease=distype;
  return(tf_mutc)
}

makefreqdf <- function(mydf,distype) {
  t_mutc=unique(mydf[,c("Hugo_Symbol","mutation","Tumor_Sample_Barcode")]);
  tf_genc=count(t_mutc, c('Hugo_Symbol')); tf_mutc=count( t_mutc, c('mutation','Hugo_Symbol')); 
  tf_mutc$gene_count<-tf_genc[match(tf_mutc$Hugo_Symbol,tf_genc$Hugo_Symbol),"freq"]; 
  tf_mutc$mut_fcount<-as.factor(tf_mutc$freq)
  tf_mutc<-tf_mutc[with(tf_mutc, order(-(gene_count),-(freq),Hugo_Symbol)),];tf_mutc$Hugo_Symbol<-factor(tf_mutc$Hugo_Symbol,levels=tf_mutc$Hugo_Symbol)
  tf_mutc$disease=distype;
  return(tf_mutc)
}