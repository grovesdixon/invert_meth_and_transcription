
library(readxl)
source('invert_meth_ge_functions.R')

read_sra_excel = function(bp){
  fn = paste('bioprojects', bp,'SraRunTable.xlsx', sep='/')
  read_excel(fn) %>% 
    select(BioProject, Run, Assay_Type)
}
sraList = map(bioprojectDirs, function(x) read_sra_excel(x))

grandSra = sraList %>% 
  purrr::reduce(rbind) 

grandSra %>%
  write_tsv(path='metadata/fullSraRunTable.txt',
            col_names = FALSE)


