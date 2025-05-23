library(mclust)
library(clustvarsel)
library(SelvarMix)
library(tidyverse)
library(Rmixmod)
library(glasso)
library(sf)

county_ag_data <- read.csv("data/county_ag_data.csv")
county_ag_scaled <- county_ag_data %>%
  ungroup() %>%
  mutate(across(!c(county, year),
                ~ (.x - mean(.x)) / (max(.x) - min(.x))))

# variable selection
# note: cut 3 variables, resulted in 5 clusters instead of 4, and a little less stable thru time
var_sel <- clustvarsel(data = county_ag_scaled[,-c(1:2)],
                       search = "greedy",
                       direction = "backward")
# another variable selection
# cuts 5 variables, results in 4 clusters
# integrates with Rmixmod package - research
var_sel2 <- SelvarClustLasso(x = county_ag_scaled[,-c(1:2)],
                             nbcluster = 1:9,
                             lambda = seq(20, 100, by=10),
                             rho = 1:5,
                             criterion = "BIC",
                             models = mixmodGaussianModel(family = "all"),
                             rmodel = c("LI", "LB", "LC"),
                             imodel = c("LI", "LB"))


#county_ag_sc_sel <- county_ag_scaled[,c("county", "year", names(var_sel$subset))]
county_ag_sc_sel <- county_ag_scaled[,c(1:2, (var_sel2$S+2))]

ag.Mclust <- Mclust(county_ag_sc_sel[,-c(1:2)])
plot(ag.Mclust, what="BIC")

ag.mclustICL <- mclustICL(county_ag_sc_sel[,-c(1:2)])
plot(ag.mclustICL)

ag_pca <- prcomp(county_ag_sc_sel[,-c(1:2)])

# factoextra::fviz_pca_biplot(ag_pca, habillage=as.factor(ag.Mclust$classification))
# factoextra::fviz_cluster(ag.Mclust)
factoextra::fviz_pca_biplot(ag_pca, habillage=as.factor(var_sel2$partition))
factoextra::fviz_cluster(var_sel2)

# county_year_ag_clust <- left_join(county_ag_data, aoi[,c("NAME", "geometry")], by = c("county" = "NAME")) %>%
#   st_as_sf() %>%
#   cbind(., ag.Mclust$classification)
# 
# ggplot(county_year_ag_clust, aes(fill = factor(ag.Mclust.classification))) +
#   geom_sf() +
#   facet_wrap(vars(year))
county_year_ag_clust <- left_join(county_ag_data, aoi[,c("NAME", "geometry")], by = c("county" = "NAME")) %>%
  st_as_sf() %>%
  cbind(., var_sel2$partition)
ggplot(county_year_ag_clust, aes(fill = factor(var_sel2.partition))) +
  geom_sf() +
  facet_wrap(vars(year))

# ag_clusts <- county_year_ag_clust %>%
#   mutate(across(-c(county, year, ag.Mclust.classification, geometry),
#                 ~scale(.x))) %>%
#   pivot_longer(cols = !c(county, year, ag.Mclust.classification, geometry),
#                names_to = "input_data_name",
#                values_to = "input_data_value")
# whisker_lengths <- ag_clusts %>%
#   st_drop_geometry() %>%
#   group_by(ag.Mclust.classification, input_data_name) %>%
#   summarise(upper_whisker = min(max(input_data_value), quantile(input_data_value, 0.75) + 1.5 * IQR(input_data_value)),
#             lower_whisker = max(min(input_data_value), quantile(input_data_value, 0.25) - 1.5 * IQR(input_data_value)),
#             mean = mean(input_data_value),
#             sd = sd(input_data_value)) %>%
#   mutate(touches_zero = ifelse((upper_whisker>0 & lower_whisker>0) | (upper_whisker<0 & lower_whisker<0), FALSE, TRUE),
#          mean_sd_zero = ifelse((mean+sd>0 & mean-sd>0) | (mean+sd<0 & mean-sd<0), FALSE, TRUE)) %>%
#   select(-c(upper_whisker, lower_whisker, mean, sd))
# 
# # actual data
# left_join(ag_clusts, whisker_lengths) %>%
#   ggplot() +
#   geom_boxplot(aes(x = input_data_value, y = input_data_name,
#                    fill = mean_sd_zero),
#                outliers = FALSE) +
#   facet_wrap(facets = "ag.Mclust.classification")
ag_clusts <- county_year_ag_clust %>%
  mutate(across(-c(county, year, var_sel2.partition, geometry),
                ~scale(.x))) %>%
  pivot_longer(cols = !c(county, year, var_sel2.partition, geometry),
               names_to = "input_data_name",
               values_to = "input_data_value")
whisker_lengths <- ag_clusts %>%
  st_drop_geometry() %>%
  group_by(var_sel2.partition, input_data_name) %>%
  summarise(upper_whisker = min(max(input_data_value), quantile(input_data_value, 0.75) + 1.5 * IQR(input_data_value)),
            lower_whisker = max(min(input_data_value), quantile(input_data_value, 0.25) - 1.5 * IQR(input_data_value)),
            mean = mean(input_data_value),
            sd = sd(input_data_value)) %>%
  mutate(touches_zero = ifelse((upper_whisker>0 & lower_whisker>0) | (upper_whisker<0 & lower_whisker<0), FALSE, TRUE),
         mean_sd_zero = ifelse((mean+sd>0 & mean-sd>0) | (mean+sd<0 & mean-sd<0), FALSE, TRUE)) %>%
  select(-c(upper_whisker, lower_whisker, mean, sd))

# actual data
left_join(ag_clusts, whisker_lengths) %>%
  ggplot() +
  geom_boxplot(aes(x = input_data_value, y = input_data_name,
                   fill = mean_sd_zero),
               outliers = FALSE) +
  facet_wrap(facets = "var_sel2.partition")

# model output
# stdev for each cluster
st_dev_clust <- as.data.frame(apply(ag.Mclust$parameters$variance$sigma, 
                                    MARGIN = 3,
                                    FUN = diag)) %>%
  mutate(across(.cols = V1:V4, .fns = sqrt)) %>%
  mutate(input_data_name = rownames(.)) %>%
  pivot_longer(cols = !input_data_name,
               names_to = "cluster",
               values_to = "stdev")

# join with mean
as.data.frame(ag.Mclust$parameters$mean) %>%
  mutate(input_data_name = rownames(.)) %>%
  pivot_longer(cols = !input_data_name,
               names_to = "cluster",
               values_to = "mean") %>%
  left_join(st_dev_clust) %>%
  filter(., (mean+stdev>0 & mean-stdev>0) | (mean+stdev<0 & mean-stdev<0)) %>%
  ggplot(aes(x = mean, y = input_data_name)) +
    geom_point() +
    geom_errorbar(aes(xmin = mean-stdev, xmax = mean+stdev)) +
    facet_wrap(vars(cluster))

barcodes <- county_year_ag_clust %>%
  select(county, year, var_sel2.partition) %>%
  pivot_wider(id_cols = c(county, geometry),
              names_from = year,
              values_from = var_sel2.partition) %>%
  rowwise() %>%
  mutate(barcode = paste(c_across(`2002`:`2022`), collapse = ""))

table(barcodes$barcode)

ggplot(barcodes, aes(fill=barcode)) +
  geom_sf()
