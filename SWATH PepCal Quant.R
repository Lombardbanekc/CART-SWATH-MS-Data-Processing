library(tidyverse)
library(fs)
library(here)
library(janitor)

##----------Protein results---------------
data <- read_delim(file.choose(), delim = '\t',
                 col_names = c('protein', 'peptide', 'mz', 'z', 'rt', 'tech1', 'tech2', 'tech3'),
                 skip = 1) %>%
  mutate(conc = 50) %>%
  rowwise() %>%
  mutate(mean = mean(c(tech1, tech2, tech3), na.rm = T),
         sd = sd(c(tech1, tech2, tech3), na.rm = T),
         rsd = sd*100/mean,
         sig_to_noise = mean/sqrt(sd)) %>%
  filter(protein != "[ RT-Cal protein ]" )

write.csv(data, file = './output/50nM_pepcal.csv')


pc_data <- tibble(filename = basename(dir_ls(path='.', glob = "*Peptides.txt"))) %>%
  mutate(file_content = map(filename, ~read_delim(.,  delim = '\t', na = c('', 'NA'),
                                                  col_names = c('protein', 'peptide', 'mz', 'z', 'rt',
                                                                'tech1', 'tech2', 'tech3'),
                                                  skip = 1))) %>%
  unnest(c(file_content))%>%
  separate('filename', c('date', 'other'), sep = '-') %>%
  separate('date', c('date', 'pc', 'conc'), sep = '_') %>%
  select(-date, -pc, -other) %>%
  filter(protein != "[ RT-Cal protein ]" ) %>%
  rowwise() %>%
  mutate(mean = mean(c(tech1, tech2, tech3), na.rm = T),
         sd = sd(c(tech1, tech2, tech3), na.rm = T),
         rsd = sd*100/mean,
         sig_to_noise = mean/sqrt(sd))

ggplot(pc_data) +
  geom_bar(aes(x = conc, y = mean), stat = 'identity') +
  facet_wrap('mz')

write_pept_results <- pc_data %>%
  group_by(mz) %>% 
  nest() %>%
  pmap(~write_csv(x = .y, path = paste0(.x, ".csv")))

##----------Peptides results----------------##

swath_pept_stat <- dir_ls(path='./output', glob = "*pepcal.csv") %>%
  map(read_csv) %>%
  reduce(bind_rows) %>%
  mutate(mean = as.numeric(mean))

write_pept_results <- swath_pept_stat %>%
  group_by(peptide) %>% 
  nest() %>%
  pmap(~write_csv(x = .y, path = paste0(.x, ".csv")))



ggplot(swath_pept_stat) +
  geom_point(aes(x = conc, y = mean)) +
facet_wrap('peptide')



