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
organisassay<-assaytx[,c(1:46,103:115)]

# veriyi normalize edelim 
organisassay_N<-log2(organisassay)
hist(organisassay_N)

# organisassay'teki sutunleri satır yapalım
transposed_data<-t(organisassay)

# veriyi PCA yapalım
dataPCA<-prcomp(transposed_data)
PC1<-dataPCA$x[,1]
PC2<-dataPCA$x[ ,2]

plot(dataPCA$x)


# ggplot(data=dataPCA$x, aes(x=PC1, y=PC2)) + geom_point()
# t test yapalım 
# once p degerlerini toplayacagimiz bir boş vektor oluşturalim

p_val=NULL


for (i in 1:nrow(organisassay_N)) { 
  p_val[i] = t.test(organisassay_N[i,1:46], organisassay_N[i,47:59])$p.value #Since t.test function gives more numbers then just the p.value, we should take the p-value by $p.value
} #Constructed a for loop for every row of the assay data. 

