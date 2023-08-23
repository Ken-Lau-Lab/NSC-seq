###Processing CSV for bulk DNA
####Processing CSV files#############################################################################################################################
################Need this file in the folder#############################################################################################################
## load libraries
library(readr)
library(Biostrings)
library(dplyr)
library(data.table)
library(reshape2)
library(stringr)

barcodes<- read.table("INUSE_PB7AllBarcodesMasterTable.txt", header = T)
barcodes<-as.data.frame(barcodes)
barcodes$spacer_new<- paste0(barcodes$Spacer, sep="", 'GGGTTAGAGCTAGAAACC')

######################################################################################################################################

args <- commandArgs(trailingOnly = TRUE)

data<- read.csv(args, header = T, na.strings=c("","NA"))


#data<- read.csv(file = '7208-MMI-1_4.csv', header = T, na.strings=c("","NA"))
data<- data[ ,c(4,14)]
colnames(data)<-c('hBC','ID')
data<- data[!is.na(data$hBC), ] 
data<- data[!is.na(data$ID), ]
data<-as.data.frame(data)

##Length sort
data$Character<- str_length(data$hBC)
data<- data[data$Character>18, ]
data$Character_Cell_ID<- str_length(data$ID)
data<- data[data$Character_Cell_ID<11, ] ##Less than 21 bp is untrimmed reads
data<- data[data$Character_Cell_ID>9, ] ##More than 12 bp is untrimmed reads

####Assign barcodes and ID
match1<- match(data$ID, barcodes$Identifier)
data$BC_Number<- barcodes$Number[match1]
data$BC<- barcodes$spacer_new[match1]
data<- data[!is.na(data$BC), ]


##Calculating ID frequency (both WT and mutated barcodes). Calculating freq before removing WT hBC.
Total_BC_Freq<- table(data$ID)
Total_BC_Freq<-as.data.frame(Total_BC_Freq)
match2<- match(data$ID, Total_BC_Freq$Var1)
data$Total_BC_Freq<- Total_BC_Freq$Freq[match2]

##Removing WT barcodes from data
data<- data[!data$hBC %in% barcodes$spacer_new, ]

#####Calculate ID frewuency for only mutated BC
Total_Mut_BC_Freq<- table(data$ID)
Total_Mut_BC_Freq<-as.data.frame(Total_Mut_BC_Freq)
match3<- match(data$ID, Total_Mut_BC_Freq$Var1)
data$Total_Mut_BC_Freq<- Total_Mut_BC_Freq$Freq[match3]

###Calculate Mutated frequency for each mutation
data$Unique<- paste0(data$hBC, sep="_", data$ID)
Mut_BC_Freq<- table(data$Unique)
Mut_BC_Freq<-as.data.frame(Mut_BC_Freq)
match4<- match(data$Unique, Mut_BC_Freq$Var1)
data$Mut_BC_Freq<- Mut_BC_Freq$Freq[match4]
data<- data[!duplicated(data[,9]), ]
data$BC_Mut_Perc<- (data$Mut_BC_Freq/data$Total_BC_Freq)*100
data<-data[,c(5,2,7,8,10,11,1,6)]

#######Mutation calling################################################################################################################################

data_join <- data %>% 
  mutate(id = seq(1,nrow(data)))

##Run the alignment
pa1 <- pairwiseAlignment(pattern = data_join$hBC , subject = data_join$BC)

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

#file = gsub(pattern = "\\.txt$", "", args)
#output_file = paste(file,'_alignment_result.txt',sep='')
## save the mutation result file
#write_csv(final_df, output_file)

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
  mutate(final_mut_pattern = paste('BC', BC_Number,'P',PatternStart,Type_new,'-',PatternSubstring_new,sep=''))



##Calculating final_mutation_pattern frequency.
final_mutation_pattern_freq<- table(data$final_mut_pattern)
final_mutation_pattern_freq<-as.data.frame(final_mutation_pattern_freq)
match5<- match(data$final_mut_pattern, final_mutation_pattern_freq$Var1)
data$final_mutation_pattern_freq<- final_mutation_pattern_freq$Freq[match5]

file = gsub(pattern = "\\.csv$", "",basename(args))
output_file = paste('output/result/',file,'_mutation.txt',sep='')
write.table(data, file = output_file, sep = "\t", row.names = F,col.names = T)


#write.table(data, file = '7208_MMI1_4_mutation.txt', sep = "\t", row.names = F,col.names = T)


