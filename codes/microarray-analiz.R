# veri okumak için paket yukleyelim 
library(GEOquery)

# verileri yukleyelim 
txdata<-getGEO("GSE50161")

# Listenin içerisindeki ilk verisetini inceleyelim 
txdata<-txdata[[1]]

#
assaytx<-txdata@assayData$exprs
phenotx<-txdata@phenoData@data
metatx<-txdata@featureData@data

#GSE50161 verisinden ependimoma ve kontrol gruplarini ayrilalim
organisassay<-assaytx[c(1:46,103:115),]

# veriyi normalize edelim 
organisassay_N<-log2(organisassay)
hist(organisassay_N)

# organisassay'teki sutunleri satır yapalım
transposed_data<-t(organisassay)

# histogram grafigi oluşturalım
hist(transposed_data)

# veriyi PCA yapalım
dataPCA<-prcomp(transposed_data)
PC1<-dataPCA$x[,1]
PC2<-dataPCA$x[ ,2]

plot(dataPCA$x)

# t test yapalım 
# once p degerlerini toplayacagimiz bir boş vektor oluşturalim

p.value <- rep(0, nrow(transposed_data))

for (i in 1:nrow(transposed_data)){
  normal<- transposed_data[i,phenotx$source_name_ch1=="normal brain"]
  tumor<- transposed_data[i, phenotx$source_name_ch1=="brain tumor"]
  
  test <- t.test(x = kontrol, y =tumor )
  p.value[i] <- test$p.value
}


sum(p.value < 0.01, na.rm = TRUE)


sum(p.value < 0.05, na.rm = TRUE)

