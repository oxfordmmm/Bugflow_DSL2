#!/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
group_size <- strtoi(args[2])
low_depth <- strtoi(args[3])
cutoff <- strtoi(args[4])
outfile <- args[5]

df <- read_tsv(file, col_names = c('id', 'locus', 'depth'))
df

ids <- df %>% 
  select(id) %>% 
  distinct() %>% 
  pull(id)

compact_df <- df %>% 
  group_by(id) %>% 
  mutate(group = ceiling(row_number() / group_size)) %>% 
  group_by(id, group) %>% 
  summarise(locus = mean(group) * group_size, depth = mean(depth)) %>% 
  ungroup() %>% 
  mutate(depth = pmin(depth, cutoff))

low_df <- compact_df %>% 
  filter(depth < low_depth)

p <- compact_df %>%  
  ggplot(aes(x=locus, y=depth)) +
  geom_point(data=low_df, color='red', size=1) +
  geom_line(linewidth=1) +
  facet_wrap(vars(id), scales='free')

ggsave(outfile, plot=p, width=25, height=12, dpi=300, limitsize = F)
