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

################################
# SES clusters

# county data
aoi <- tigris::counties("ID") %>%
  select(c("GEOID", "NAME"))

ru_data <- read.csv("data/RU_data.csv")
rs_data <- read.csv("data/RS_data.csv")
resource <- left_join(ru_data, rs_data)

ses_data_list <- list(Resouce = resource)

ses_data_list_sc <- lapply(ses_data_list, function(df) {
  df %>%
    na.omit() %>%
    ungroup() %>%
    mutate(across(!c(county, year),
                  ~ (.x - mean(.x)) / (max(.x) - min(.x))))
})

ses_var_sel <- lapply(ses_data_list_sc, function(df) {
  SelvarClustLasso(x = df[,-c(1:2)],
                   nbcluster = 1:9,
                   lambda = seq(20, 100, by=10),
                   rho = 1:5,
                   criterion = "BIC",
                   models = mixmodGaussianModel(family = "all"),
                   rmodel = c("LI", "LB", "LC"),
                   imodel = c("LI", "LB"))
})

ses_data_sel <- mapply(function(df, varsel) {
  df[,c(1:2, (varsel$S+2))]
}, df=ses_data_list_sc, varsel=ses_var_sel)

# ses_Mclust <- lapply(ses_data_sel, function(df) {
#   Mclust(df[,-c(1:2)], G=1:15)
# })
# 
# ses_mclustICL <- lapply(ses_data_sel, function(df) {
#   mclustICL(df[,-c(1:2)], G=1:15)
# })

# visualize cluster structure and principal components
pca_viz <- lapply(ses_data_sel, function(df) {
  prcomp(df[,-c(1:2)])
})
# p1 <- factoextra::fviz_pca_biplot(pca, habillage=as.factor(mod$classification))
# p2 <- factoextra::fviz_cluster(mod)
factoextra::fviz_pca_biplot(pca_viz[[1]], habillage=as.factor(ses_var_sel[[1]]$partition),
                            title = "Resource Units")
factoextra::fviz_pca_biplot(pca_viz[[2]], habillage=as.factor(ses_var_sel[[2]]$partition),
                            title = "Resource System")

# map of clusters in all counties, all years
clust_dat <- mapply(function(df, mod, system_name) {
  left_join(df, aoi[,c("NAME", "geometry")], by = c("county" = "NAME")) %>%
    st_as_sf() %>%
    cbind(., mod$partition) %>%
    #mutate(cluster = paste0("R", mod.partition)) %>%
    select(-mod.partition)
}, df = ses_data_list_sc, mod = ses_var_sel, system_name=names(ses_data_list_sc))
ggplot(clust_dat[["Resource"]], aes(fill = cluster)) +
  geom_sf() +
  facet_wrap(vars(year))

# visualize data difference between clusters
clust_dat_pv <- lapply(clust_dat, function(df) {
  df %>%
    pivot_longer(cols = !c(county, year, cluster, geometry),
                 names_to = "input_data_name",
                 values_to = "input_data_value")
})
whiskers <- lapply(clust_dat_pv, function(df) {
  df %>%
    st_drop_geometry() %>%
    group_by(cluster, input_data_name) %>%
    summarise(upper_whisker = min(max(input_data_value), quantile(input_data_value, 0.75) + 1.5 * IQR(input_data_value)),
              lower_whisker = max(min(input_data_value), quantile(input_data_value, 0.25) - 1.5 * IQR(input_data_value)),
              mean = mean(input_data_value),
              sd = sd(input_data_value)) %>%
    mutate(touches_zero = ifelse((upper_whisker>0 & lower_whisker>0) | (upper_whisker<0 & lower_whisker<0), FALSE, TRUE),
           mean_sd_zero = ifelse((mean+sd>0 & mean-sd>0) | (mean+sd<0 & mean-sd<0), FALSE, TRUE)) %>%
    select(-c(upper_whisker, lower_whisker, mean, sd))
})
# viz with actual data
left_join(clust_dat_pv[["RS"]], whiskers[["RS"]]) %>%
  ggplot() +
  geom_boxplot(aes(x = input_data_value, y = input_data_name,
                   fill = mean_sd_zero),
               outliers = FALSE) +
  facet_wrap(facets = "cluster")

barcodes <- bind_rows(clust_dat) %>%
  select(county, year, cluster) %>%
  group_by(county, geometry, year) %>%
  mutate(yearcode = paste0(cluster, collapse = "")) %>%
  select(-cluster) %>%
  unique() %>%
  pivot_wider(id_cols = c(county, geometry),
              names_from = year,
              values_from = yearcode) %>%
  select(county, `1990`, `2000`, `2010`, `2020`, geometry) %>%
  rowwise() %>%
  mutate(barcode = paste(c_across(`1990`:`2020`), collapse = "-"))

# almost no common barcodes
table(barcodes$barcode)

ggplot(barcodes, aes(fill=barcode)) +
  geom_sf()

system_barcodes <- lapply(clust_dat, function(df) {
  df %>%
    select(county, year, cluster) %>%
    pivot_wider(id_cols = c(county, geometry),
                names_from = year,
                values_from = cluster) %>%
    select(county, `1990`, `2000`, `2010`, `2020`, geometry) %>%
    rowwise() %>%
    mutate(barcode = paste(c_across(`1990`:`2020`), collapse = "-"))
})
table(system_barcodes[["RS"]]$barcode)
ggplot(system_barcodes[["RU"]], aes(fill=barcode)) +
  geom_sf()
