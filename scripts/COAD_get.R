library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg38)

#data source
#https://api.gdc.cancer.gov/data/d28f95fc-af3b-497d-806b-eb0625d6c831

#COAD
tissueName = 'COAD'
aav1898_Data_S2_tab1 <- read.csv("aav1898_Data_S2_tab1.csv", stringsAsFactors=FALSE)
tissueSpecificDF = aav1898_Data_S2_tab1[grepl(tissueName,aav1898_Data_S2_tab1$Peak_Name),]
colnames(tissueSpecificDF)[1] = 'Chromosome'

fileConn<-file("COAD_data.fa",'w')
for(i in 1:nrow(aav1898_Data_S2_tab1)){
  writeLines(paste(">",aav1898_Data_S2_tab1$Peak_Name[i]), fileConn)
  tempSeq = getSeq(Hsapiens, aav1898_Data_S2_tab1$Chromosome[i], start=aav1898_Data_S2_tab1$Start[i], end=aav1898_Data_S2_tab1$End[i], strand="+", as.character=TRUE)
  writeLines(tempSeq, fileConn)
}
close(fileConn)

