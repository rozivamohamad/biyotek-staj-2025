library(GEOquery)

txdata<-getGEO("GSE12345")


exprs_data <- exprs(txdata)
 # Ekspresyon matrisi
pdata <- pData(txdata) 
 # Fenotipik bilgiler (örneğin grup, hastalık durumu)
head(txdata[,1:5])
# İlk 5 örneğin ilk genleri
head(pdata)