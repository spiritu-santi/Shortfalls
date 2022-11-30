library(KnowBR)
library(tidyverse)
load("output/pams_for_completeness.Rdata")
bd <- for_comp$tenth_degree %>% filter(LEVEL2_NAM %in% c("Mexico", "Central America")) %>%
  select(c("Species", "Longitude",  "Latitude", "Counts")) %>% 
  #slice(1:10000) %>% 
  #filter(Longitude%in%c(-113.5,-112.5,-110.5,-109.5,-103.5,-102.5,-101.5,-84.5),Latitude%in%c(29.5)) %>% 
  as.data.frame()

data=bd


bd %<>% as_tibble() %>% unite("cell",Longitude:Latitude,sep="_") %>% group_by(Species,cell) %>% summarize(Counts=sum(Counts)) %>% 
  pivot_wider(names_from = Species, values_from = Counts) %>% select(-cell)

bd

cu <- vegan::specaccum(bd, method = "rarefaction")

cumodelo <- try(nls(cu.richness ~ (A + B * cu.sites)/(1 + C * cu.sites), data = datosc, trace = T, start = list(A = 1, B = 1, C = 0)), silent = TRUE)

data("adworld")
KnowB(data=bd)
