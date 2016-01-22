suppressWarnings(suppressMessages(library(argparser,logical.return=FALSE,verbose=FALSE,quietly=TRUE)))
suppressWarnings(suppressMessages(library(reshape2,logical.return=FALSE,verbose=FALSE,quietly=TRUE)))
suppressWarnings(suppressMessages(library(ggplot2,logical.return=FALSE,verbose=FALSE,quietly=TRUE)))
suppressWarnings(suppressMessages(library(Matrix,logical.return=FALSE,verbose=FALSE,quietly=TRUE)))

mycols=c("n_ref_count","n_alt_count","t_ref_count","t_alt_count")
col.pos <- c(1,2,3,5,6,7,9,10,14,16,17,match(mycols, scan(input_consolidated,sep="\t",what=character(0),nlines=1)))

indir="/home/jgrewal/projects/tcga_lowvaf/augmaf_results/consolidated/"
hnsc=paste(indir,"HNSC.txt",sep="");
luad=paste(indir,"LUAD.txt",sep=""); 
lusc=paste(indir,"LUSC_trimmed.txt",sep=""); 

readmaf <- function(filepath,distype,deflist=c(1,2,3,5,6,7,9,10,14,16,17)){
  mycols=c("n_ref_count","n_alt_count","t_ref_count","t_alt_count")
  col.pos <- c(deflist,match(mycols, scan(filepath,sep="\t",what=character(0),nlines=1)))
  df <- ((read.table(pipe( paste( "cut -d$'\t' -f",paste(col.pos,collapse=',')," ",filepath,sep='')), header=TRUE, sep = "\t",as.is=TRUE)))
  df["nvaf"]=(df["n_alt_count"])/(df["n_ref_count"]+df["n_alt_count"])
  df["tvaf"]=(df["t_alt_count"])/(df["t_ref_count"]+df["t_alt_count"])
  df$mutation = paste(df$Hugo_Symbol,df$Chromosome, df$Start_Position, df$End_Position, sep="_")
  df$disease = distype; #df=df[!is.na(df$nvaf),]
  return(df)
}
hnsc_df=readmaf(hnsc,"HNSC")
lusc_df=readmaf(lusc,"LUSC",deflist=c(1:11)); lusc_df$dbSNP_RS=ifelse(grepl('^rs',lusc_df$dbSNP_RS),lusc_df$dbSNP_RS,"novel")
luad_df=readmaf(luad,"LUAD"); luad_df$dbSNP_RS=ifelse(grepl('^rs',luad_df$dbSNP_RS),luad_df$dbSNP_RS,"novel")
df<- na.omit(rbind(hnsc_df,luad_df,lusc_df))

#Unrestricted data
clean_table=df[(df$nvaf == 0),!names(df) %in% c("n_ref_count","n_alt_count","t_ref_count","t_alt_count")]
clean_matrix=clean_table[,c("Hugo_Symbol","tvaf","Tumor_Sample_Barcode")]
clean_matrix=acast(clean_matrix,Hugo_Symbol~Tumor_Sample_Barcode,value.var="tvaf",fun.aggregate = mean)
nnzero(clean_matrix,na.counted=TRUE)
clean_matrix[is.na(clean_matrix)] <- 0 #Replace all na's with 0
binary_matrix <- clean_matrix; binary_matrix[binary_matrix > 0] <- 1

#Biclustering by presence of 1's, biclust
suppressWarnings(suppressMessages(library(biclust,logical.return=FALSE,verbose=FALSE,quietly=TRUE)))
c_value=nnzero(clean_matrix,na.counted=TRUE)/(dim(clean_matrix)[2] * dim(clean_matrix)[1])
rCclust <- biclust(binary_matrix, method=BCCC(),delta=c_value,alpha=NA,number=7 )
rplaid <- biclust(binary_matrix, method=BCPlaid(),row.release=0.5, col.release=0.5, max.layers=7 )
#Uses SVD. Try with "irrc" idependent rescaling of rows and columns
rSpectral <- biclust(clean_matrix,method=BCSpectral(), normalization="bistochastization",withinVar=0.1) 

#Biclustering sparse binary genomic data (BicBin)



#Some interesting variant classes
indir=paste(indir,"results/",sep="")
write.table(df[(df$tvaf < df$nvaf)&(df$tvaf>0),],file=paste(indir,"lowt_highn_genes.txt",sep=""),append=FALSE,quote=FALSE,row.names=FALSE,sep="\t")

#Generate low freq novel variants' table
choosenovel <- function(df,threshold=0.5,stringent=0){
  df=df[df$Variant_Type=="SNP"  & df$tvaf>0 & (df$nvaf==0),] ; df_highvaf=df[df$tvaf>=threshold,]
  df_lowvaf=df[df$tvaf<threshold & !(df$mutation %in% df_highvaf$mutation),] #mutation filter needed? & !(df$Hugo_Symbol %in% df_highvaf$Hugo_Symbol)
  if(stringent){
  df_lowvaf = df[df$tvaf<threshold & !(df$Hugo_Symbol %in% df_highvaf$Hugo_Symbol),]}
  df_novel_lowvaf=df_lowvaf[df_lowvaf$dbSNP_RS=="novel",]
  #df_novel=df[df$dbSNP_RS=="novel",]
  return(df_novel_lowvaf)
}
hnsc_novel_lowvaf=choosenovel(hnsc_df);
lowvaf=0.5;
write.table(choosenovel(df,threshold=lowvaf), file=paste(indir,"novel_",lowvaf,"vaf_mutations.txt",sep=""),append=FALSE,quote=FALSE,row.names=FALSE,sep="\t")
write.table(choosenovel(df,threshold=lowvaf,stringent=1), file=paste(indir,"novel_",lowvaf,"vaf_genes.txt",sep=""),append=FALSE,quote=FALSE,row.names=FALSE,sep="\t")
mut_slow=choosenovel(df,threshold=lowvaf);gene_slow=choosenovel(df,threshold=lowvaf,stringent=1);

#Get counts
library(plyr);library(reshape2); library(gplots); library(RColorBrewer)
pdf(paste(indir,"gene_lowvaf_hist.pdf"),height=35,width=20)
gene_hmap = count(gene_slow,c("Hugo_Symbol","disease")); gene_hmap$dfreq=gene_hmap$freq/length(unique(gene_slow$Tumor_Sample_Barcode)); 
gene_hmapd=acast(gene_hmap[gene_hmap$freq>1,],Hugo_Symbol~disease,value.var='freq'); gene_hmapd[is.na(gene_hmapd)]<-0; 
heatmap.2(na.omit(gene_hmapd),cellnote=gene_hmapd, lwid=c(1.5,5),lhei=c(1.5,14),tracecol=NA,dendrogram="row",cexCol=1,cexRow=1,margins=c(4,16),
          main="Pan cohort Low vaf genes, patient frequency across cohorts", col=colorRampPalette(c("white","green","blue2","violet","purple"))(100), notecol="black")
dev.off()
write.table(gene_hmapd,file=paste(indir,"gene_lowvaf_hist_data.txt",sep=""),append=FALSE,quote=FALSE,row.names=TRUE,sep="\t")

mut_hmap = count(mut_slow,c("mutation","disease")); mut_hmap$dfreq=mut_hmap$freq/length(unique(mut_slow$Tumor_Sample_Barcode)); 
mut_hmapd=acast(mut_hmap[mut_hmap$freq>2,],mutation~disease,value.var='freq');mut_hmapd[is.na(mut_hmapd)]<-0; 
heatmap.2((mut_hmapd),tracecol=NA,dendrogram="row",cexCol=1,cexRow=1,margins=c(5,16))


#Mutation frequency matrix
makefreqdf <- function(mydf,distype,mutcols=c('Hugo_Symbol','Tumor_Sample_Barcode')) {
  t_mutc=unique(mydf[,c("Hugo_Symbol","mutation","Tumor_Sample_Barcode")]);
  tf_genc=count(t_mutc, c('Hugo_Symbol')); tf_mutc=count( t_mutc, mutcols); 
  tf_mutc$gene_count<-tf_genc[match(tf_mutc$Hugo_Symbol,tf_genc$Hugo_Symbol),"freq"]; 
  tf_mutc$mut_fcount<-as.factor(tf_mutc$freq)
  tf_mutc<-tf_mutc[with(tf_mutc, order(-(gene_count),-(freq),Hugo_Symbol)),];tf_mutc$Hugo_Symbol<-factor(tf_mutc$Hugo_Symbol,levels=tf_mutc$Hugo_Symbol)
  tf_mutc$disease=distype;
  return(tf_mutc)
}

# Low freq genes
'''{r}
library(reshape2); library(plyr); library(ggplot2)
hnsc_pdf = makefreqdf(choosenovel(hnsc_df,lowvaf,1),"HNSC")
luad_pdf=makefreqdf(choosenovel(luad_df,lowvaf,1),"LUAD")
lusc_pdf=makefreqdf(choosenovel(lusc_df,lowvaf,1),"LUSC")
tf_mutc<- na.omit(rbind(hnsc_pdf,luad_pdf,lusc_pdf))

ggplot(tf_mutc,aes(x=Hugo_Symbol,y=freq,fill=gene_count)) + geom_bar(stat='identity',position="stack",subset=.(freq>0),) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=1));
ggplot(tf_mutc[tf_mutc$freq>1 & tf_mutc$gene_count>1,],aes(x=Hugo_Symbol,fill=freq)) + geom_bar(stat='bin',position="stack") +  theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=1)) + geom_text(aes(x=Hugo_Symbol, y=-1,label=gene_count),stat='identity',size=2) +
  scale_y_continuous(breaks=seq(0,25),1) + facet_grid(disease ~ .)

ggplot(tf_mutc[(tf_mutc$Hugo_Symbol %in% tf_mutc[(tf_mutc$gene_count>1),"Hugo_Symbol"]),],aes(x=Hugo_Symbol,y=freq,fill=mut_fcount)) + geom_bar(stat='identity',position="stack",subset=.(freq>40)) + theme(text = element_text(size=30),axis.text.x=element_text(angle=90,hjust=1,vjust=1)) +
  labs(title="Cohort specific Low variant genes, > 1 patient") + facet_grid(disease ~ ., scales="free_y") +
  geom_text(aes(x=Hugo_Symbol, y=-1,label=gene_count),stat='identity',size=8)+ guides(fill = guide_legend(keywidth = 3, keyheight = 6));
write.table(tf_mutc[(tf_mutc$Hugo_Symbol %in% tf_mutc[(tf_mutc$gene_count>1),"Hugo_Symbol"]),],file=paste(indir,"summary-gene-lowvaf-",lowvaf,"-cohort_data.txt",sep=""),append=FALSE,quote=FALSE,row.names=FALSE,sep="\t")
ggsave(file=paste(indir,"/summary-gene-lowvaf-",lowvaf,"-cohort.pdf",sep=""),width=200,height=40,limitsize=FALSE)
'''
#Pan cohort
hnsc_pdf = makefreqdf(gene_slow[gene_slow$disease=="HNSC",],"HNSC") 
luad_pdf=makefreqdf(gene_slow[gene_slow$disease=="LUAD",],"LUAD")
lusc_pdf=makefreqdf(gene_slow[gene_slow$disease=="LUSC",],"LUSC")
tf_mutc2<- na.omit(rbind(hnsc_pdf,luad_pdf,lusc_pdf))

ggplot(tf_mutc2,aes(x=Hugo_Symbol,y=freq,fill=gene_count)) + geom_bar(stat='identity',position="stack",subset=.(freq>0),) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=1));
ggplot(tf_mutc2[tf_mutc2$freq>1 & tf_mutc2$gene_count>1,],aes(x=Hugo_Symbol,fill=gene_count)) + geom_bar(stat='bin',position="stack") +  theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=1)) + geom_text(aes(x=Hugo_Symbol, y=-1,label=gene_count),stat='identity',size=2) +
  scale_y_continuous(breaks=seq(0,25),1) + facet_grid(disease ~ .)

pancoh_genlow=tf_mutc2[(tf_mutc2$Hugo_Symbol %in% tf_mutc2[(tf_mutc2$gene_count>1),"Hugo_Symbol"]),]
ggplot(pancoh_genlow,aes(x=Hugo_Symbol,y=freq,fill=mut_fcount)) + geom_bar(stat='identity',position="stack",subset=.(freq>0)) + theme(text = element_text(size=30),axis.text.x=element_text(angle=90,hjust=1,vjust=1)) +
  labs(title="Pan cohort Low variant genes, > 1 patient, pan cohort filter") + facet_grid(disease ~ ., scales="free_y") +
  geom_text(aes(x=Hugo_Symbol, y=-1,label=gene_count),stat='identity',size=8)+ guides(fill = guide_legend(keywidth = 3, keyheight = 6));
ggsave(file=paste(indir,"/summary-gene-lowvaf-all-pancohort.pdf",sep=""),width=49,height=40)
write.table(pancoh_genlow,file=paste(indir,"summary-gene-lowvaf-pancohort_data.txt",sep=""),append=FALSE,quote=FALSE,row.names=FALSE,sep="\t")

mydat=tf_mutc[(tf_mutc$Hugo_Symbol %in% tf_mutc[(tf_mutc$gene_count>1),"Hugo_Symbol"]),]
mydat=mydat[!(mydat$Hugo_Symbol %in% pancoh_genlow$Hugo_Symbol), ]
ggplot(mydat,aes(x=Hugo_Symbol,y=freq,fill=mut_fcount)) + geom_bar(stat='identity',position="stack",subset=.(freq>0)) + theme(text = element_text(size=30),axis.text.x=element_text(angle=90,hjust=1,vjust=1)) +
  labs(title="Cohort specific Low variant genes, > 2 patients \n excludes consistently low freq genes") + facet_grid(disease ~ ., scales="free_y") +
  geom_text(aes(x=Hugo_Symbol, y=-1,label=gene_count),stat='identity',size=8)+ guides(fill = guide_legend(keywidth = 3, keyheight = 6));
write.table(mydat,file=paste(indir,"filtered-gene-lowvaf-cohort_atleast3_data.txt",sep=""),append=FALSE,quote=FALSE,row.names=FALSE,sep="\t")
ggsave(file=paste(indir,"/filtered-gene-lowvaf-cohort_atleast3.pdf",sep=""),width=49,height=40,limitsize=FALSE)

#Mutation uniqueness based - low freq mutations only

library(reshape2); library(plyr); library(ggplot2)
hnsc_pdf = makefreqdf(choosenovel(hnsc_df,lowvaf,0),"HNSC",mutcols=c('mutation','Hugo_Symbol'))
luad_pdf=makefreqdf(choosenovel(luad_df,lowvaf,0),"LUAD",mutcols=c('mutation','Hugo_Symbol'))
lusc_pdf=makefreqdf(choosenovel(lusc_df,lowvaf,0),"LUSC",mutcols=c('mutation','Hugo_Symbol'))
tf_mutc<- na.omit(rbind(hnsc_pdf,luad_pdf,lusc_pdf))
tf_aggr<- aggregate(tf_mutc[,c("Hugo_Symbol","freq","gene_count","disease")],by=list("Hugo_Symbol","disease"),FUN=mean,na.rm=TRUE)
#df_aggr <- with(hnsc_pdf, aggregate (hnsc_pdf, by = list(Hugo_Symbol), FUN = mean))
ggplot(tf_mutc[tf_mutc$gene_count>5,],aes(x=Hugo_Symbol,y=freq,fill=mut_fcount)) + geom_bar(stat='identity',position="stack",subset=.(freq>0),) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=1));
ggplot(tf_mutc[tf_mutc$freq>1 & tf_mutc$gene_count>1,],aes(x=Hugo_Symbol,fill=mut_fcount)) + geom_bar(stat='bin',position="stack") +  theme_bw() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=1)) + geom_text(aes(x=Hugo_Symbol, y=-1,label=gene_count),stat='identity',size=2) +
    scale_y_continuous(breaks=seq(0,25),1) + facet_grid(disease ~ .)

ggplot(tf_mutc[(tf_mutc$Hugo_Symbol %in% tf_mutc[(tf_mutc$freq>2),"Hugo_Symbol"]) & (tf_mutc$Hugo_Symbol == "TTN"),],aes(x=mutation,y=freq,fill=mut_fcount)) + geom_bar(stat='identity',position="stack",subset=.(freq>0)) + theme(text = element_text(size=30),axis.text.x=element_text(angle=90,hjust=1,vjust=1)) +
  labs(title="Cohort specific Low variant mutations, > 2 patients \n excludes consistently low freq genes") + facet_grid(disease ~ ., scales="free_y") +
  geom_text(aes(x=Hugo_Symbol, y=-1,label=gene_count),stat='identity',size=8)+ guides(fill = guide_legend(keywidth = 3, keyheight = 6));
ggsave(file=paste(indir,"ttn-mutations-lowvaf-atleast3-cohort.pdf",sep=""),width=49,height=20,limitsize=FALSE)

ggplot(tf_mutc[tf_mutc$gene_count>25,],aes(x=Hugo_Symbol,y=gene_count,fill=disease)) + geom_bar(stat='identity',subset=.(freq>1)) + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=1)) +
  labs(title="Low variant genes, "); 
outdir="/Users/jgrewal/Desktop/"
ggsave(file=paste(indir,"/summary-3-novel-lowvaf-hnsc.pdf",sep=""),width=49,height=40)


#DEPRECATED code

ggplot(df_lowvaf,aes(x=Tumor_Sample_Barcode)) + geom_histogram() + labs(title="SNVs per sample", x = "Tumour Sample", y= "SNV count")
ggplot(df_novel,aes(x=Tumor_Sample_Barcode)) + geom_histogram() + labs(title="Novel SNVs per sample", x = "Tumour Sample", y= "SNV count")
df2=df_lowvaf
ggplot(df_lowvaf,aes(x=Hugo_Symbol)) + geom_bar(stat="bin")+ labs(title="Mutation frequency per gene", x = "Gene", y= "SNV count")
ggplot(data=df2, aes(x=tvaf, fill=Tumor_Sample_Barcode)) + geom_density(alpha=0.5) + labs(title="Density distribution of SNVs per patient", y = "SNV density", x= "VAF_tumour")
ggplot(data=df2, aes(x=tvaf, colour=Tumor_Sample_Barcode)) + geom_freqpoly() + labs(title="Frequency distribution of SNVs per patient", y = "SNV count", x= "VAF_tumour")


#Plot high freq gene counts
gfreq=data.frame(table(df_novel_lowvaf$Hugo_Symbol))
ggplot(gfreq[gfreq$Freq>20,],aes(x=Var1,y=Freq)) + geom_point() + theme(axis.text.x=element_text(angle=90,hjust=1,vjust=1)) + facet_grid(. ~ disease)
outdir="/Users/jgrewal/Desktop/"
ggsave(file=paste(outdir,"/genefreq-novel-hnsc.pdf",sep=""),width=35,height=20)
