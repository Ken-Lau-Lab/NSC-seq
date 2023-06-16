## load library
library(readr)
library(Biostrings)
library(dplyr)
library(data.table)
library(reshape2)

#args <- commandArgs(trailingOnly = TRUE)
data <- read_delim('7984-MMI-1_1_R2_trimmed_mutation.txt', 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE)
data<- read.csv("7118_MMI1_Mito.csv", header = T, na.strings=c("","NA"))
data<-data[!is.na(data$reads), ]
data<- data[!is.na(data$BC), ]


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
  mutate(final_mut_pattern = paste('BC',BC,'P',PatternStart,Type_new,'-',PatternSubstring_new,sep=''))

### for single cell only
data_matrix<-dcast(data, data$final_mut_pattern ~ data$Cell, fun.aggregate = function(x){as.integer(length(x) > 0)})
data_matrix<-as.data.frame(data_matrix)
write.table(data, file = "7118_MMI1_MC1_mitoBC_mutation.txt", sep = "\t", row.names = T,col.names = T)
