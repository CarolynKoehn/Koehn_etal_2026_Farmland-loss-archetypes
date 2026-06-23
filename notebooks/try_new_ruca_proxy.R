library(tidyverse)

dist_from_urb <- read.csv("~/Documents/Chapter_2/FarmPrES-archetypes/data/clean_data_for_PCA/urb_dist_areas.csv")
commutes <- read.csv("~/Documents/Chapter_2/FarmPrES-archetypes/data/clean_data_for_PCA/commutes_all.csv")

commutes2 <- commutes %>%
  pivot_longer(names_to = "cat",
               cols = !c(county, year)) %>%
  mutate(commute = word(cat, 2, sep="_")) %>%
  select(-cat) %>%
  group_by(county, year, commute) %>%
  summarise(commute_perc = sum(value)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(county,year),
              names_from = "commute",
              values_from = "commute_perc")

urbs_plus_commutes <- full_join(dist_from_urb, commutes2, by=c("county", "year"))

upc_pca <- prcomp(urbs_plus_commutes[,-c(1:2)])
