library(tidyverse)
library(fs)
library(here)
library(janitor)

##----------Protein results----------------##

swath_prot_files <- files <- dir_ls(path = '.', glob = '*Proteins.txt') 


swath_prot_stat <- tibble(filename = basename(dir_ls(path='.', glob = "*Proteins.txt"))) %>%
  mutate(file_content = map(filename, ~read_delim(.,  delim = '\t', na = c('', 'NA'), 
                                                  col_names = c('protein', 'tech1', 'tech2', 'tech3'),
                                                  skip = 1))) %>%
  unnest(c(file_content)) %>%
  clean_names() %>%
    rowwise() %>%
  mutate(mean = mean(c(tech1, tech2, tech3), na.rm = T),
         log_mean = log(mean),
         sd = sd(c(tech1, tech2, tech3), na.rm = T),
         rsd = sd*100/mean) %>%
  filter(protein != "[ RT-Cal protein ]")


write_prot_results <- swath_prot_stat %>%
  group_by(filename) %>% 
  nest() %>%
  pmap(~write_csv(x = .y, path = paste0(.x, "_processed.csv")))

n_protein <- swath_prot_stat %>%
  group_by(filename) %>%
  summarise(n = n()) %>%
  write_csv('n_protein_summary.csv')

n_protein_graph <- swath_prot_stat %>%
  group_by(filename) %>%
  summarise(n = n()) 

ggplot(n_protein_graph) +
  geom_bar(aes(x = filename, y = n), stat = 'identity')

n_20cv_protein <- swath_prot_stat %>%
  filter(rsd<20) %>%
  group_by(filename) %>%
  summarise(n = n()) %>%
  write_csv('n_protein_20cv_summary.csv')

ggplot(n_20cv_protein) +
  geom_bar(aes(x = filename, y = n), stat = 'identity')

##----------Peptides results----------------##

swath_pept_stat <- tibble(filename = basename(dir_ls(path='.', glob = "*Peptides.txt"))) %>%
  mutate(file_content = map(filename, ~read_delim(.,  delim = '\t', na = c('', 'NA'),
                                                  col_names = c('protein', 'peptide', 'mz', 'z', 'rt',
                                                                'tech1', 'tech2', 'tech3'),
                                                  skip = 1))) %>%
  unnest(c(file_content)) %>%
  clean_names() %>%
  rowwise() %>%
  mutate(mean = mean(c(tech1, tech2, tech3), na.rm = T),
         log_mean = log(mean),
         sd = sd(c(tech1, tech2, tech3), na.rm = T),
         rsd = sd*100/mean) %>%
  filter(protein != "[ RT-Cal protein ]")

write_pept_results <- swath_pept_stat %>%
  group_by(filename) %>% 
  nest() %>%
  pmap(~write_csv(x = .y, path = paste0(.x, "_processed.csv")))

n_peptides <- swath_pept_stat %>%
  group_by(filename) %>%
  summarise(n = n()) %>%
  write_csv('n_peptide_summary.csv')

n_20cv_peptides <- swath_pept_stat %>%
  filter(rsd<20) %>%
  group_by(filename) %>%
  summarise(n = n()) %>%
  write_csv('n_20percCV_peptide_summary.csv')

ggplot(n_peptides) +
  geom_bar(aes(x = filename, y = n), stat = 'identity')



