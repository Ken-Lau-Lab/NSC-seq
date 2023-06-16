## load library
library(readr)
library(Biostrings)
library(dplyr)
library(data.table)
library(reshape2)
library(qdapRegex)
library(stringr)
#####################################Processing CSV files and aggregating###############################

BC_mito <- read.table(file = 'INUSE_BC_Mito_12.txt', header = T, fill = F)
#BC_mito<- BC_mito[which(BC_mito$Gene=='Atp6' | BC_mito$Gene=='Co2'| BC_mito$Gene=='Cyb'), ] ###No need to run this if there are 12 mito genes

mito_gene = BC_mito$Gene

for (gene in mito_gene){
  BC_mito_filter = BC_mito[BC_mito$Gene==gene,]
  filenames = Sys.glob(file.path('output_mito', paste0("*mito_",gene,"\\.csv")))
  data <- read.csv(filenames)
  if (nrow(data)>0){
    data<- data %>% 
    select(c(11,21)) %>% 
    filter_at(vars(X0_R2, Cell), all_vars(!is.na(.))) %>% 
    mutate(reads = unlist(rm_between(X0_R2, BC_mito_filter$Start, BC_mito_filter$End, extract=TRUE))) %>% 
    mutate(Character = str_length(reads)) %>% 
    mutate(Character_Cell_ID = str_length(Cell)) %>% 
    mutate(BC = BC_mito_filter$BC) %>% 
    filter(!is.na(reads))
  
  median_length<- mean(data$Character, na.rm=T)
  median_length<- round(median_length, digits = 0)
  
  data<- data[!data$Character> (median_length + 2), ] 
  data<- data[!data$Character< (median_length - 2), ]
  
  data<- data[data$Character_Cell_ID<21, ] ##Less than 21 bp is untrimmed reads
  data<- data[data$Character_Cell_ID>12, ] ##More than 12 bp is untrimmed reads
  
  data_mut<- data[!data$reads %in% BC_mito_filter$BC, ]
  if (nrow(data_mut)> 0){
    data_mut$Gene<- BC_mito_filter$Gene
    data_mut$Cell<- paste0(data_mut$Cell, sep="-", 'test-1')
    output<- data_mut[, c(3,2,6,7)]
    write_csv(output,paste0("MMI1/",gene,'_test-1.csv'))
  }
  }
  
}


filenames = Sys.glob(file.path('MMI1','*.csv'))

MMI1_Mito <- filenames %>% 
  lapply(read.csv) %>% 
  bind_rows


MMI1_Mito$Unique<- paste0(MMI1_Mito$reads, sep="_", MMI1_Mito$Cell)

Unique<- table(MMI1_Mito$Unique)
Unique<-as.data.frame(Unique)
match_mito<- match(MMI1_Mito$Unique, Unique$Var1)
MMI1_Mito$Freq<- Unique$Freq[match_mito]
MMI1_Mito<- MMI1_Mito[!duplicated(MMI1_Mito[,5]), ]

#write_csv(MMI1_Mito,'8337_MMI2_Mito.csv')

########################################################################Mutation Calling####################################


#args <- commandArgs(trailingOnly = TRUE)
#data1 <- read_csv(file = '8337_MMI6_Mito.csv')
#data2 <- read_csv(file = '8337_MMI12_Mito.csv')
#data<-rbind(data1,data2)

data<- MMI1_Mito

#data <- read.table(file = 'data_7995_MMI1.txt', header = T, fill = F)
#data<-data[!is.na(data$V1), ]
#data<- data[!is.na(data$BC), ]


data_join <- data %>% 
  mutate(id = seq(1,nrow(data)))

##Run the alignment
pa1 <- pairwiseAlignment(pattern = data_join$reads , subject = data_join$BC)

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
  mutate(final_mut_pattern = paste(Gene,'P',PatternStart,Type_new,'-',PatternSubstring_new,sep=''))



write.table(data, file = "Test_mitoBC_table.txt", sep = "\t", row.names = F,col.names = T)

