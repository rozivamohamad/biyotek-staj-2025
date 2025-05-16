# veri okumak için paket yukleyelim.
library(GEOquery)
library(ggplot2)

# verileri yukleyelim 
txdata<-getGEO("GSE50161")

# Listenin içerisindeki ilk verisetini inceleyelim 
txdata<-txdata[[1]]

#
assaytx<-txdata@assayData$exprs
phenotx<-txdata@phenoData@data
metatx<-txdata@featureData@data

#GSE50161 verisinden ependimoma ve kontrol gruplarini ayrilalim
organisassay<-assaytx[,c(1:46,103:115)]

#histogram grafigi olusturalim
hist(organisassay)



# veriyi normalize edelim 
organisassay_N <-log2(organisassay)
hist(organisassay_N)

# organisassay'teki sutunleri satır yapalım
transposed_data<-t(organisassay)



# veriyi PCA yapalım
dataPCA<-prcomp(transposed_data)
PC1<-dataPCA$x[,1]
PC2<-dataPCA$x[ ,2]



labels<-row.names(transposed_data)
ggplot(data = dataPCA$x, aes(x= PC1, y= PC2)) + geom_point() + geom_text(label=labels, check_overlap = T)

# t test yapalım 
# once p degerlerini toplayacagimiz bir boş vektor oluşturalim

p_value <- NULL

for (i in 1:nrow(organisassay_N)){

  p_value[i]<-t.test(organisassay_N [i,1:46], organisassay_N [i,47:59])$p.value
  
  
}

BF_p_value<-p.adjust(p_value,method = "bonferroni")

BF_sign_genes=which(BF_p_value<0.01)

#gen isimleri geri dondurelim.
genes_names<-metatx$`Gene Symbol`[BF_sign_genes]

unique_genisimler<-unique(genes_names)

##dondurdumuz gen isinleri bir txt yapalım.
write.table(unique_genisimler,file = "unique gen isimleri.txt",sep="\t",quote = F,row.names = F)


BF_p_val<- p.adjust(p_val, method = "bonferroni")
BF_sign_genes = which(BF_p_val < 0.01)

Sign_Gene_Names = metatx$`Gene Symbol`[BF_sign_genes]
Unique_Gene_Names = unique(Sign_Gene_Names)
write.table(Unique_Gene_Names, file = "Unique Gene Names.txt", sep = "\t", quote = F, row.names = F)

# Sample Clustering

Clust_Samples = organisassay_N[BF_sign_genes,]
dist_samples = as.dist(1-cor(Clust_Samples, method = "spearman")) #Calculating the distance between samples using pearson method
hSamples = hclust(dist_samples, method = "complete") #Hierarchically clustering the genes using complete linkage method.
plot(as.dendrogram(hSamples)) #plotting the dendrogram"

#Gene Clustering

Clust_Genes = t(organisassay_N[BF_sign_genes,]) #Taking the transpose of the data, since gene clustering is required.
dist_genes = as.dist(1-cor(Clust_Genes, method = "spearman")) #Finding the inter cluster distance between genes.
hGenes = hclust(dist_genes, method = "complete") #Clustering the samples using above distance and complete linkage method.
plot(as.dendrogram(hGenes)) #Plotting the clustered genes as dendrogram.

heatmap(Clust_Samples, Rowv = as.dendrogram(hGenes), Colv = as.dendrogram(hSamples), scale = "row",) #Creating a heatmap
