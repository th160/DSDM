library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)
setwd("C:\\Users\\thchen7\\Downloads\\coding")

#data source
#https://api.gdc.cancer.gov/data/d28f95fc-af3b-497d-806b-eb0625d6c831

#LUSC
tissueName = 'LUSC'
aav1898_Data_S2_tab1 <- read.csv("aav1898_Data_S2_tab1.csv", stringsAsFactors=FALSE)
tissueSpecificDF = aav1898_Data_S2_tab1[grepl(tissueName,aav1898_Data_S2_tab1$Peak_Name),]
colnames(tissueSpecificDF)[1] = 'Chromosome'


fileConn<-file("LUSC_data.fa",'w')
for(i in 1:nrow(tissueSpecificDF)){
  writeLines(paste(">",tissueSpecificDF$Peak_Name[i]), fileConn)
  tempSeq = getSeq(Hsapiens, tissueSpecificDF$Chromosome[i], start=tissueSpecificDF$Start[i], end=tissueSpecificDF$End[i], strand="+", as.character=TRUE)
  writeLines(tempSeq, fileConn)
}
close(fileConn)
