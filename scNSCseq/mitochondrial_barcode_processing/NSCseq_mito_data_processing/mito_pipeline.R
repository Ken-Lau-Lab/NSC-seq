library(qdapRegex)
library(readr)
library(stringr)
library(dplyr)

BC_mito <- read.table(file = 'INUSE_BC_Mito_12.txt', header = T, fill = F)

mito_gene = BC_mito$Gene

for (gene in mito_gene){
  BC_mito_filter = BC_mito[BC_mito$Gene==gene,]
  filenames = Sys.glob(file.path('output_mito2', paste0("*-[1]*",gene,"\\.csv")))
  data <- filenames %>% 
    lapply(read.csv) %>% 
    bind_rows %>% 
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
  if (nrow(data_mut)!=0){
    data_mut$Gene<- BC_mito_filter$Gene
    data_mut$Cell<- paste0(data_mut$Cell, sep="-", '7118-MMI1')
    
    output<- data_mut[, c(3,2,6,7)]
    write_csv(output,paste0("MMI2/",gene,'_7118_MMI1.csv'))
  }
}


filenames = Sys.glob(file.path('MMI2','*.csv'))

MMI1_Mito <- filenames %>% 
  lapply(read.csv) %>% 
  bind_rows


MMI1_Mito$Unique<- paste0(MMI1_Mito$reads, sep="_", MMI1_Mito$Cell)

Unique<- table(MMI1_Mito$Unique)
Unique<-as.data.frame(Unique)
match_mito<- match(MMI1_Mito$Unique, Unique$Var1)
MMI1_Mito$Freq<- Unique$Freq[match_mito]
MMI1_Mito<- MMI1_Mito[!duplicated(MMI1_Mito[,5]), ]

write_csv(MMI1_Mito,'7118_MMI1_Mito.csv')
