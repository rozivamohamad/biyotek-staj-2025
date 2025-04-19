library(GEOquery)

txdata<-getGEO("GSE48350")


txdata<-txdata[[1]]

assaytx<-txdata@assayData$exprs

phenotx<-txdata@phenoData@data

metatx<-txdata@featureData@data

transposed_data<-t(assaytx)
hist(transposed_data)

transposed_data_N<-log2(transposed_data)

hist(transposed_data_N)
