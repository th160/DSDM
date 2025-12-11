library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(biomaRt)
library(DNAshapeR)
library(depmixS4)

#data source
#https://genome-euro.ucsc.edu/cgi-bin/hgTables?db=hg19&hgta_group=regulation&hgta_track=wgEncodeRegDnaseClustered&hgta_table=wgEncodeRegDnaseClusteredV3&hgta_doSchema=describe+table+schema

centerLen = 100
wgEncodeRegDnaseClusteredV3 <- read.delim("D:/Coding/GRF2020/wgEncodeRegDnaseClusteredV3.bed", header=FALSE)
#AllCellTypemask = (wgEncodeRegDnaseClusteredV3$V4>=125)

filePath = "Dnase-seq_data.fa"
fileConn<-file(filePath,'w')
for(i in 1:nrow(wgEncodeRegDnaseClusteredV3)){
  if(wgEncodeRegDnaseClusteredV3[i,4]>=125){
    seqName = paste(wgEncodeRegDnaseClusteredV3[i,1],wgEncodeRegDnaseClusteredV3[i,2],wgEncodeRegDnaseClusteredV3[i,3],sep='-')
    writeLines(paste(">",seqName), fileConn)
    midpoint = round((wgEncodeRegDnaseClusteredV3[i,2]+wgEncodeRegDnaseClusteredV3[i,3])/2)
    tempSeq = getSeq(Hsapiens, wgEncodeRegDnaseClusteredV3[i,1], start=midpoint-(centerLen/2), end=midpoint+(centerLen/2)-1, strand="+", as.character=TRUE)
    writeLines(tempSeq, fileConn)
  } 
}
close(fileConn)

