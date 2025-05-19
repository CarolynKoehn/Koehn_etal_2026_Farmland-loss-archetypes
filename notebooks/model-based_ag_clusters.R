library(mclust)

ag.Mclust <- Mclust(county_ag_scaled[,-c(1:2)])
plot(ag.Mclust, what="BIC")
# chose 4 clusters with VVE

ag.mclustICL <- mclustICL(county_ag_scaled[,-c(1:2)])
plot(ag.mclustICL)
# chose 4 clusters with VVE

ag_pca <- prcomp(county_ag_scaled[,-c(1:2)])

factoextra::fviz_pca_biplot(ag_pca, habillage=as.factor(ag.Mclust$classification))
factoextra::fviz_cluster(ag.Mclust)

county_year_ag_clust <- left_join(county_ag_data, aoi[,c("NAME", "geometry")], by = c("county" = "NAME")) %>%
  st_as_sf() %>%
  cbind(., ag.Mclust$classification)

ggplot(county_year_ag_clust, aes(fill = factor(ag.Mclust.classification))) +
  geom_sf() +
  facet_wrap(vars(year))

ag_clusts <- county_year_ag_clust %>%
  mutate(across(-c(county, year, ag.Mclust.classification, geometry),
                ~scale(.x))) %>%
  pivot_longer(cols = !c(county, year, ag.Mclust.classification, geometry),
               names_to = "input_data_name",
               values_to = "input_data_value")
whisker_lengths <- ag_clusts %>%
  st_drop_geometry() %>%
  group_by(ag.Mclust.classification, input_data_name) %>%
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
  facet_wrap(facets = "ag.Mclust.classification")

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
  select(county, year, ag.Mclust.classification) %>%
  pivot_wider(id_cols = c(county, geometry),
              names_from = year,
              values_from = ag.Mclust.classification) %>%
  rowwise() %>%
  mutate(barcode = paste(c_across(`2002`:`2022`), collapse = ""))

table(barcodes$barcode)

ggplot(barcodes, aes(fill=barcode)) +
  geom_sf()
