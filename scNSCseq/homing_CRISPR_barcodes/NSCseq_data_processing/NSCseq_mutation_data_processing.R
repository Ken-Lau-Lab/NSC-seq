#####Processing CSV files#############################################################################################################################
##Set output folder as current directory in R

#data1<- read.csv("8019-MMI-3.csv", header = T, na.strings=c("","NA"))
#data2<- read.csv("8019-MMI-6.csv", header = T,na.strings=c("","NA"))
#data<-rbind(data1,data2)
data<- read.csv("output/Test.csv", header = T, na.strings=c("","NA"))
data<- data[ ,c(12,21)]
data<- data[!is.na(data$hBC_R2), ] 
data<- data[!is.na(data$Cell), ]
data<-as.data.frame(data)
##Length sort
library("stringr")
data$Character<- str_length(data$hBC_R2)
data<- data[data$Character>25, ]
data$Character_Cell_ID<- str_length(data$Cell)
data<- data[data$Character_Cell_ID<21, ] ##Less than 21 bp is untrimmed reads
data<- data[data$Character_Cell_ID>12, ] ##More than 12 bp is untrimmed reads


################Need this file in the folder#############################################################################################################
barcodes<- read.table("INUSE_PB7AllBarcodesMasterTable.txt", header = T)
barcodes<-as.data.frame(barcodes)
barcodes$spacer_new<- paste0(barcodes$Spacer, sep="", 'GGGTTAGAGCTAGAA')

data<- data[!data$hBC_R2 %in% barcodes$spacer_new, ]
data$Unique<- paste0(data$hBC_R2, sep="_", data$Cell)

Freq<- table(data$Unique)
Freq<-as.data.frame(Freq)
match1<- match(data$Unique, Freq$Var1)
data$Freq<- Freq$Freq[match1]
data<- data[!duplicated(data[,5]), ]
data<-data[,c(1,2,6)]

##Adding project name with cell ID as it would be like in scRNAseq count matrix###########################################################################################################
data$Cell<- paste0(data$Cell, sep="-", 'test-1')

#########################
data$BC<- NA
data_MMI = data.frame()
for (i in seq(9,6)){
  data<- data[is.na(data$BC), ]
  barcodes$bp<- str_sub(barcodes$Spacer, 1, i)
  data$bp<- str_sub(data$hBC_R2, 1, i)
  match2<- match(data$bp, barcodes$bp)
  data$BC<- barcodes$Number[match2]
  data_match<- data[!is.na(data$BC), ]
  data_MMI = rbind(data_MMI,data_match)
}
#5
data<- data[is.na(data$BC), ]
barcodes$bp<- str_sub(barcodes$Spacer, 1, 5)
data$bp<- str_sub(data$hBC_R2, 1, 5)
match2<- match(data$bp, barcodes$bp)
data$BC<- barcodes$Number[match2]
data_bp5<- data[!is.na(data$BC), ]
data_bp5<-data_bp5[data_bp5$BC== c(10,28,34), ]

data_MMI<- rbind(data_MMI,data_bp5)

match3<- match(data_MMI$BC, barcodes$Number)
data_MMI$spacer_new<- barcodes$spacer_new[match3]

#write.table(data_MMI, file = "data_test.txt", sep = "\t", row.names = F,col.names = T)

#######Mutation calling################################################################################################################################
## load libraries
library(readr)
library(Biostrings)
library(dplyr)
library(data.table)
library(reshape2)

data<- data_MMI
data<- data[!is.na(data$BC), ]

data_join <- data %>% 
  mutate(id = seq(1,nrow(data)))

##Run the alignment
pa1 <- pairwiseAlignment(pattern = data_join$hBC_R2 , subject = data_join$spacer_new)

pa2 <- pairwiseAlignment(pattern = as.character(aligned(pattern(pa1)))
                         , subject = as.character(aligned(subject(pa1)))
)

result <- mismatchTable(pa2)

result <- result %>% select(PatternId,PatternStart,PatternSubstring,SubjectSubstring) %>% 
  mutate(Type = case_when(
    PatternSubstring == '-' ~ 'deletion',
    SubjectSubstring == '-' ~ 'insertion',
    TRUE ~ 'mismatch')
  ) %>% 
  left_join(data_join,by = c('PatternId' = 'id'))

table <- result

new_df = data.frame()
for (id in unique(table$PatternId)) {
  temp <- table %>% filter(PatternId == id)
  
  old_start <- NA
  starts = temp$PatternStart[1]
  start_index = 1
  inorder = TRUE
  for (start in temp$PatternStart) {
    if(!is.na(old_start)) {
      if(start != old_start+1 | temp[temp$PatternStart==start,]$Type != temp[temp$PatternStart==old_start,]$Type) { # & type is the same TODO
        starts[[length(starts)+1]] <- start
        start_index[[length(start_index)+1]] <- which(start == temp$PatternStart)
      }
    }
    old_start = start
  }
  for(i in 1:length(start_index)) {
    if(i != length(start_index)) mutation <- temp$PatternStart[start_index[i]]:temp$PatternStart[start_index[i+1]-1]
    else {
      mutation <- temp$PatternStart[start_index[i]]:tail(temp$PatternStart,n=1)
    }
    
    subset_df <- temp[temp$PatternStart %in% mutation,]
    
    new_row <- subset_df[1,]
    if(new_row$Type == "insertion") new_row$PatternSubstring = paste0(subset_df$PatternSubstring, collapse='')
    if(new_row$Type == "deletion") new_row$SubjectSubstring = paste0(subset_df$SubjectSubstring, collapse='')
    if(new_row$Type == "mismatch") {
      new_row$SubjectSubstring = paste0(subset_df$SubjectSubstring, collapse='')
      new_row$PatternSubstring = paste0(subset_df$PatternSubstring, collapse='')
    }
    
    if(nrow(new_df) == 0) new_df = new_row
    else new_df <- rbind(new_df, new_row)
  }
}

## sort the result
final_df <- new_df[with(new_df, order(PatternId, PatternStart)),]


data <- final_df
data<- data %>% mutate(Type_new = case_when(
  Type=='insertion'~ "I",
  Type=='deletion'~"D",
  Type=='mismatch'~"M"
))%>%
  mutate(PatternSubstring_new = case_when(
    PatternSubstring=='-' ~ SubjectSubstring,
    TRUE ~ PatternSubstring
  )) %>%
  mutate(final_mut_pattern = paste('BC', BC,'P',PatternStart,Type_new,'-',PatternSubstring_new,sep=''))


write.table(data, file = "Test_mutations_table.txt", sep = "\t", row.names = F,col.names = T)

