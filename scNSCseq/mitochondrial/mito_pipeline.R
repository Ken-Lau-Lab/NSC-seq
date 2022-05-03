library(qdapRegex)
library(readr)
BC_mito <- read_delim("INUSE_BC_Mito.txt", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)

mito_gene = BC_mito$Gene

for (gene in mito_gene){
  BC_mito_filter = BC_mito[BC_mito$Gene==gene,]
  ## extract all -3 and -6 files
  filenames = Sys.glob(file.path('output_mito', paste0("*-[3,6]*",gene,"\\.csv")))
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
  data_mut$Gene<- BC_mito_filter$Gene
  data_mut$Cell<- paste0(data_mut$Cell, sep="-", '7995-MMI3')
  
  output<- data_mut[, c(3,2,6,7)]
  ## create a new directory MMI3 to save all output files
  write_csv(output,paste0("MMI3/",gene,'_7995_MMI3.csv'))
  
}


filenames = Sys.glob(file.path('MMI3','*.csv'))

MMI3_Mito <- filenames %>% 
  lapply(read.csv) %>% 
  bind_rows


MMI3_Mito$Unique<- paste0(MMI3_Mito$reads, sep="_", MMI3_Mito$Cell)

Unique<- table(MMI3_Mito$Unique)
Unique<-as.data.frame(Unique)
match_mito<- match(MMI3_Mito$Unique, Unique$Var1)
MMI3_Mito$Freq<- Unique$Freq[match_mito]
MMI3_Mito<- MMI3_Mito[!duplicated(MMI3_Mito[,5]), ]

write_csv(MMI3_Mito,'MMI3_Mito.csv')
