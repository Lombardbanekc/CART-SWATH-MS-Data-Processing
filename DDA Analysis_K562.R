library(tidyverse)
library(fs)
library(openxlsx)

#Set Working Directory
File_path <- setwd("C:/Users/lombardbanekc/Documents/09 DATA/03 Processed Data/2020 Columns Evaluation/K562 IDA_new")

Data_xl <- read.xlsx(file.choose(), sheet = 'Search Summary', startRow = 4)


## Protein results ##
prot_0.01FDR <- Data_xl %>%
                filter( Data.Level == 'Protein', FDR.Type == 'Global', FDR == 0.01)

Prot_Data <-  read.xlsx(file.choose(), sheet = 'Protein Summary')

Prot_Res<- Prot_Data %>% 
    mutate(rownames = as.numeric(row.names(.)), expt = 'Ek-new_1') %>%
    filter(rownames <= prot_0.01FDR[1,"ID.Yield"], .$'Peptides(95%)' >= 1) %>%
    select(N, Total, '%Cov', '%Cov(50)', '%Cov(95)', Accession, Name, `Peptides(95%)`, expt) %>%
    write.csv('K562_IDA_Proteins_Ek-new_1.csv')


## Peptides results ##
pept_0.01FDR <- Data_xl %>%
    filter( Data.Level == 'Distinct peptide', FDR.Type == 'Global', FDR == 0.01)

Pept_Data <-  read.xlsx(file.choose(), sheet = 'Distinct Peptide Summary')
Pept_Res<- Pept_Data %>% 
    mutate(rownames = as.numeric(row.names(.)), expt = 'PX_COMB') %>%
    filter(rownames <= pept_0.01FDR[1,"ID.Yield"]) %>%
    write.csv('K562_IDA_Peptides_PX_COMB.csv')
