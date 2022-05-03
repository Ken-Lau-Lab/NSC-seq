library(stringr)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)

barcodes<- read.table("INUSE_PB7AllBarcodesMasterTable.txt", header = T)
barcodes<-as.data.frame(barcodes)
barcodes$spacer_new<- paste0(barcodes$Spacer, sep="", 'GGGTTAGAGCTAGAAATAGCAA')

data<- read.table(args, header=FALSE,fill=TRUE,stringsAsFactor=FALSE)
data$Character<- str_length(data$V1)
data<- data[data$Character<71 & data$Character>24, ]
match1<- match(data$V1, barcodes$spacer_new)
data$BC<- barcodes$Number[match1]

data_MMI = data.frame()
for (i in seq(9,6)){
  data<- data[is.na(data$BC), ]
  barcodes$bp<- str_sub(barcodes$Spacer, 1, i)
  data$bp<- str_sub(data$V1, 1, i)
  match2<- match(data$bp, barcodes$bp)
  data$BC<- barcodes$Number[match2]
  data_match<- data[!is.na(data$BC), ]
  data_MMI = rbind(data_MMI,data_match)
}

#5
data<- data[is.na(data$BC), ]
barcodes$bp<- str_sub(barcodes$Spacer, 1, 5)
data$bp<- str_sub(data$V1, 1, 5)
match2<- match(data$bp, barcodes$bp)
data$BC<- barcodes$Number[match2]
data_bp5<- data[!is.na(data$BC), ]
data_bp5<-data_bp5[data_bp5$BC== c(10,28,34), ]

data_MMI<- rbind(data_MMI,data_bp5)

match3<- match(data_MMI$BC, barcodes$Number)
data_MMI$spacer_new<- barcodes$spacer_new[match3]
data_unique<- table(data_MMI$V1)
data_unique<-as.data.frame(data_unique)
match4<- match(data_MMI$V1, data_unique$Var1)
data_MMI$Freq<-data_unique$Freq[match4]
data_MMI<- data_MMI[!duplicated(data_MMI[,1]), ]

file = gsub(pattern = "\\.txt$", "", args)
output_file = paste(file,'_mutation.txt',sep='')

write.table(data_MMI, file = output_file, sep = "\t", row.names = F,col.names = T)