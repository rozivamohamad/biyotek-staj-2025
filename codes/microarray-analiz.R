# veri okumak için paket yukleyelim 
library(GEOquery)

# verileri yukleyelim 
txdata<-getGEO("GSE48350")

# Listenin içerisindeki ilk verisetini inceleyelim 
txdata<-txdata[[1]]

#
assaytx<-txdata@assayData$exprs
phenotx<-txdata@phenoData@data
metatx<-txdata@featureData@data

# assaytx'teki sutunleri satır yapalım
transposed_data<-t(assaytx)

# histogram grafigi oluşturalım
hist(transposed_data)

# veriyi normalize edelim 
transposed_data_N<-log2(transposed_data)

# histogram grafigi olusturalım 
hist(transposed_data_N)
