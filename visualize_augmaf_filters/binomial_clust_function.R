library(vegan)
shannon_measure = diversity(clean_matrix,index="simpson")

library(amap)
data.xfr = choosenovel(df,threshold=0.5,stringent=1)
clusters=Kmeans(x=data.xfr[,c("Hugo_Symbol","tvaf")],centers=20, method="pearson") #,"Tumor_Sample_Barcode"
