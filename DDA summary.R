library(tidyverse)
library(stats)
library(readr)

## Protein results ##
files <- dir_ls(path = './output', glob = '*proteins_COMB*')

  R <- files %>% map(read_csv) %>%
  reduce(bind_rows) %>%
  mutate(separate(., col = expt, into = c('expt', 'replicate'), sep = '_', remove = T)) %>%
  group_by(expt, replicate) %>%
  summarise(n_proteins = n())


ggplot(R) +
  geom_bar(aes(x = expt, y = n_proteins), stat = 'identity')
 


write.csv(R, 'output/proteins_summary_COMB.csv')
  

## Peptide results ##
  
files_pept <- dir_ls(path = './output', glob = '*peptides*')
  
pept_sum <- files_pept %>% map(read_csv) %>%
    reduce(bind_rows) %>%
    mutate(separate(., col = expt, into = c('expt', 'replicate'), sep = '_', remove = T)) %>%
    group_by(expt, replicate) %>%
    summarise(n_peptides = n()) %>%
    group_by(expt) %>%
    summarise(mean = mean(n_peptides), SD = sd(n_peptides))
  
  ggplot(pept_sum) +
    geom_bar(aes(x = expt, y = mean), stat = 'identity') +
    geom_errorbar(aes(x = expt, ymin = mean-SD, ymax = mean + SD), stat = 'identity', width = 0.3)  

  write.csv(pept_sum, 'output/peptidess_summary.csv')
  
