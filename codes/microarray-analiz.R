library(GEOquery)

txdata<-getGEO("GSE48350")


txdata<-txdata[[1]]

assaytx<-txdata@assayData$exprs

phenotx<-txdata@phenoData@data

metatx<-txdata@featureData@data

hippocampus_rows<-phenotx[grep("hippocampus",phenotx$characteristics_ch1.1),]

hippocampus_rows_AD<-phenotx[grep("hippocampus",phenotx$characteristics_ch1.2),]

postcentral_gyrus_rows<-phenotx[grep("postcentral gyrus",phenotx$characteristics_ch1.1),]

superior_frontal_gyrus_rows<-phenotx[grep("superior frontal gyrus",phenotx$characteristics_ch1.1),]

entorhinal_cortex_rows<-phenotx[grep("entorhinal cortex",phenotx$characteristics_ch1.1),]


united_hippocampus<-rbind(hippocampus_rows,hippocampus_rows_AD)

hippocampus_data <- assaytx[, colnames(assaytx) %in% united_hippocampus$geo_accession]


hist(hippocampus_data)

hippocampus_data_N<-log2(hippocampus_data)

hist(hippocampus_data_N)