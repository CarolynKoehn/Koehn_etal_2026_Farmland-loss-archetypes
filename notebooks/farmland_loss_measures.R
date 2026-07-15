library(tidyverse)
library(terra)
library(sf)
library(tigris)

conv_01_16 <- rast("R:/SarahHalperin/Data/AFT_Data/ag_land_conversion_2001_2016_ID_30m (1).tif")
lulc_01 <- rast("R:/SarahHalperin/Data/AFT_Data/land_cover_and_use_2001_ID_30m (1).tif")
lulc_16 <- rast("R:/SarahHalperin/Data/AFT_Data/land_cover_and_use_2016_ID_30m (1).tif")
pvr <- rast("R:/SarahHalperin/Data/AFT_Data/productivity_versatility_resiliency_2016_conus_120m.tif")

aoi <- tigris::counties("ID", progress_bar = FALSE) %>%
  select(c("GEOID", "NAME")) %>%
  filter(NAME %in% c("Jefferson", "Gem", "Madison", "Clark", "Jerome", "Bonneville", "Teton", "Power", "Blaine", "Elmore", "Butte", "Adams", "Minidoka", "Bingham", "Ada", "Cassia", "Payette", "Fremont", "Bannock", "Twin Falls", "Owyhee", "Washington", "Gooding", "Lincoln", "Camas", "Canyon")) %>% st_transform(st_crs(lulc_01))

mean_pvr_aoi <- mean(values(crop(pvr, st_transform(aoi, st_crs(pvr)), mask=TRUE)), na.rm=TRUE)

ada <- counties("ID") %>% filter(NAME == "Ada") %>% st_transform(st_crs(lulc_01))
ada01 <- crop(lulc_01, ada, mask=TRUE)
ada16 <- crop(lulc_16, ada, mask=TRUE)

table(paste(as.data.frame(ada01, xy=TRUE)$first,
            as.data.frame(ada16, xy=TRUE)$first,
            sep="-"))

# get counties with cluster labels
clust_dat <- list()
clust_dat[["Resource"]] <- read.csv("data/model results/clust_dat_resource.csv") %>%
  left_join(aoi, by = c("county" = "NAME")) %>%
  st_as_sf()
clust_dat[["Actors"]] <- read.csv("data/model results/clust_dat_actors.csv") %>%
  left_join(aoi, by = c("county" = "NAME")) %>%
  st_as_sf()
county_clusts <- left_join(st_drop_geometry(clust_dat[["Resource"]][,c("county", "year", "cluster")]),
                           st_drop_geometry(clust_dat[["Actors"]][,c("county", "year", "cluster")]),
                           by = c("county", "year"),
                           suffix = c(".R", ".A")) %>%
  mutate(cluster.R = factor(cluster.R,
                            levels = c("R1", "R4", "R5", "R2", "R3"),
                            labels = c("R1", "R2", "R3", "R4", "R5")),
         cluster.A = factor(cluster.A,
                            levels = c("A1", "A2", "A5", "A3", "A4"),
                            labels = c("A1", "A2", "A3", "A4", "A5")))
county_barcodes <- county_clusts %>%
  rowwise() %>%
  mutate(comb_clust = paste0(cluster.R,cluster.A)) %>%
  pivot_wider(id_cols="county",
              names_from = year,
              values_from = comb_clust)
county_barcodes$barcode <- do.call(paste, c(county_barcodes[,-1], sep="-"))
county_barcodes$barcode_A <- str_remove_all(county_barcodes$barcode, "R[:digit:]")
county_barcodes$barcode_R <- str_remove_all(county_barcodes$barcode, "A[:digit:]")

county_barcodes[,c("crop_to_UHD", "crop_to_UHD_avg_pvr",
                   "crop_to_LDR", "crop_to_LDR_avg_pvr",
                   "past_to_UHD", "past_to_UHD_avg_pvr",
                   "past_to_LDR", "past_to_LDR_avg_pvr",
                   "range_to_UHD", "range_to_UHD_avg_pvr",
                   "range_to_LDR", "range_to_LDR_avg_pvr")] <- NA
for(i in 1:nrow(county_barcodes)) {
  lulc01c <- crop(lulc_01, subset(aoi, NAME==county_barcodes$county[i]), mask=TRUE)
  lulc16c <- crop(lulc_16, subset(aoi, NAME==county_barcodes$county[i]), mask=TRUE)
  pvrc <- crop(pvr, subset(aoi, NAME==county_barcodes$county[i]), mask=TRUE)
  
  lulc_change <- lulc01c
  values(lulc_change) <- paste(as.character(values(lulc01c)), as.character(values(lulc16c)), sep="-")
  lulc_change_counts <- freq(lulc_change)
  
  # 1 pixel = 900 m2 = 0.2223948433 acres
  county_barcodes$crop_to_UHD[i] <- sum(
    lulc_change_counts[which(lulc_change_counts$value %in% c("1-6", "1-1106")),"count"]
    )*0.2223948433
  county_barcodes$crop_to_UHD_avg_pvr[i] <- 
    ifelse(test = county_barcodes$crop_to_UHD[i] == 0,
           yes = NA,
           no = mean(extract(x = pvrc,
                             y = as.points(ifel(lulc_change %in% c("1-6", "1-1106"), 1, NA)),
                             fun=mean, na.rm=TRUE)[,2],
                     na.rm=TRUE))
  
  county_barcodes$crop_to_LDR[i] <- sum(
    lulc_change_counts[which(lulc_change_counts$value %in% c("1-1101", "1-1102", "1-1103", "1-1104", "1-1105", "1-1106", "1-1107", "1-1108", "1-1109", "1-1110", "1-1112")),"count"]
    )*0.2223948433
  county_barcodes$crop_to_LDR_avg_pvr[i] <- 
    ifelse(test = county_barcodes$crop_to_LDR[i] == 0,
           yes = NA, 
           no = mean(extract(x = pvrc,
                             y = as.points(ifel(lulc_change %in% c("1-1101", "1-1102", "1-1103", "1-1104", "1-1105", "1-1106", "1-1107", "1-1108", "1-1109", "1-1110", "1-1112"), 1, NA)),
                             fun=mean, na.rm=TRUE)[,2],
                     na.rm=TRUE))

  county_barcodes$past_to_UHD[i] <- sum(
    lulc_change_counts[which(lulc_change_counts$value %in% c("2-6", "2-1106")),"count"]
  )*0.2223948433
  county_barcodes$past_to_UHD_avg_pvr[i] <- 
    ifelse(test = county_barcodes$past_to_UHD[i] == 0,
           yes = NA,
           no = mean(extract(x = pvrc,
                             y = as.points(ifel(lulc_change %in% c("2-6", "2-1106"), 1, NA)),
                             fun=mean, na.rm=TRUE)[,2],
                     na.rm=TRUE))
  
  county_barcodes$past_to_LDR[i] <- sum(
    lulc_change_counts[which(lulc_change_counts$value %in% c("2-1101", "2-1102", "2-1103", "2-1104", "2-1105", "2-1106", "2-1107", "2-1108", "2-1109", "2-1110", "2-1112")),"count"]
  )*0.2223948433
  county_barcodes$past_to_LDR_avg_pvr[i] <- 
    ifelse(test = county_barcodes$past_to_LDR[i] == 0,
           yes = NA,
           no = mean(extract(x = pvrc,
                             y = as.points(ifel(lulc_change %in% c("2-1101", "2-1102", "2-1103", "2-1104", "2-1105", "2-1106", "2-1107", "2-1108", "2-1109", "2-1110", "2-1112"), 1, NA)),
                             fun=mean, na.rm=TRUE)[,2],
                     na.rm=TRUE))
  
  county_barcodes$range_to_UHD[i] <- sum(
    lulc_change_counts[which(lulc_change_counts$value %in% c("3-6", "3-1106")),"count"]
  )*0.2223948433
  county_barcodes$range_to_UHD_avg_pvr[i] <- 
    ifelse(test = county_barcodes$range_to_UHD[i] == 0,
           yes = NA,
           no = mean(extract(x = pvrc,
                             y = as.points(ifel(lulc_change %in% c("3-6", "3-1106"), 1, NA)),
                             fun=mean, na.rm=TRUE)[,2],
                     na.rm=TRUE))
    
  county_barcodes$range_to_LDR[i] <- sum(
    lulc_change_counts[which(lulc_change_counts$value %in% c("3-1101", "3-1102", "3-1103", "3-1104", "3-1105", "3-1106", "3-1107", "3-1108", "3-1109", "3-1110", "3-1112")),"count"]
  )*0.2223948433
  county_barcodes$range_to_LDR_avg_pvr[i] <- 
    ifelse(test = county_barcodes$range_to_LDR[i] == 0,
           yes = NA,
           no = mean(extract(x = pvrc,
                             y = as.points(ifel(lulc_change %in% c("3-1101", "3-1102", "3-1103", "3-1104", "3-1105", "3-1106", "3-1107", "3-1108", "3-1109", "3-1110", "3-1112"), 1, NA)),
                             fun=mean, na.rm=TRUE)[,2],
                     na.rm=TRUE))
  
  print(i)
}

county_barcodes$barcode_A_lab <- ifelse(test = str_extract(county_barcodes$barcode_A, "^A[:digit:]") == str_extract(county_barcodes$barcode_A, "A[:digit:]$"),
                                        yes = str_extract(county_barcodes$barcode_A, "^A[:digit:]"),
                                        no = paste(str_extract(county_barcodes$barcode_A, "^A[:digit:]"), str_extract(county_barcodes$barcode_A, "A[:digit:]$"), sep="\u279E"))

county_barcodes$barcode_R_lab <- ifelse(test = str_extract(county_barcodes$barcode_R, "^R[:digit:]") == str_extract(county_barcodes$barcode_R, "R[:digit:]$"),
                                        yes = str_extract(county_barcodes$barcode_R, "^R[:digit:]"),
                                        no = paste(str_extract(county_barcodes$barcode_R, "^R[:digit:]"), str_extract(county_barcodes$barcode_R, "R[:digit:]$"), sep="\u279E"))

# exploratory plots
ggplot(county_barcodes, aes(x=crop_to_UHD, y=crop_to_UHD_avg_pvr, shape=barcode_R_lab, color = barcode_A_lab)) + geom_point()
ggplot(county_barcodes, aes(x=crop_to_LDR, y=crop_to_LDR_avg_pvr, shape=barcode_R_lab, color = barcode_A_lab)) + geom_point()
ggplot(county_barcodes, aes(x=past_to_UHD, y=past_to_UHD_avg_pvr, shape=barcode_R_lab, color = barcode_A_lab)) + geom_point()
ggplot(county_barcodes, aes(x=past_to_LDR, y=past_to_LDR_avg_pvr, shape=barcode_R_lab, color = barcode_A_lab)) + geom_point()
ggplot(county_barcodes, aes(x=range_to_UHD, y=range_to_UHD_avg_pvr, shape=barcode_R_lab, color = barcode_A_lab)) + geom_point()
ggplot(county_barcodes, aes(x=range_to_LDR, y=range_to_LDR_avg_pvr, shape=barcode_R_lab, color = barcode_A_lab)) + geom_point()

county_barcodes %>%
  rowwise() %>%
  mutate(ag_to_UHD = sum(c_across(c(crop_to_UHD, past_to_UHD, range_to_UHD)), na.rm=T),
         ag_to_LDR = sum(c_across(c(crop_to_LDR, past_to_LDR, range_to_LDR)), na.rm=T)) %>%
  ggplot(aes(x=ag_to_UHD, y=ag_to_LDR, shape=barcode_R_lab, color = barcode_A_lab)) +
  geom_abline(slope=1, intercept = 0, linetype="dashed") +
  geom_point()

county_barcodes %>%
  rowwise() %>%
  mutate(pvr_to_UHD = mean(c_across(c(crop_to_UHD_avg_pvr, past_to_UHD_avg_pvr, range_to_UHD_avg_pvr)), na.rm=T),
         pvr_to_LDR = mean(c_across(c(crop_to_LDR_avg_pvr, past_to_LDR_avg_pvr, range_to_LDR_avg_pvr)), na.rm=T)) %>%
  ggplot(aes(x=pvr_to_UHD, y=pvr_to_LDR, shape=barcode_R_lab, color = barcode_A_lab)) +
  geom_abline(slope=1, intercept = 0, linetype="dashed") +
  geom_point()

# these are the ones I think (maybe without R1-A1?)
plot_UHD_data <- county_barcodes %>%
  rowwise() %>%
  mutate(ag_to_UHD = sum(c_across(c(crop_to_UHD, past_to_UHD, range_to_UHD)), na.rm=T),
         pvr_to_UHD = mean(c_across(c(crop_to_UHD_avg_pvr, past_to_UHD_avg_pvr, range_to_UHD_avg_pvr)), na.rm=T))
  
ggplot(aes(x=ag_to_UHD, y=pvr_to_UHD, shape=barcode_R_lab, color = barcode_A_lab),
       data= plot_UHD_data) +
  geom_hline(yintercept = mean_pvr_aoi, linetype="dashed") +
  geom_text(aes(x=7100, y=0.263, label="Avg Study Area PVR"), hjust="right", inherit.aes = F, data=data.frame()) +
  geom_segment(data=subset(plot_UHD_data, county %in% c("Teton", "Canyon", "Twin Falls")), 
            inherit.aes = F,
            aes(x=ag_to_UHD+200, xend=ag_to_UHD,
                y=pvr_to_UHD-0.01, yend=pvr_to_UHD)) +
  geom_text(data=subset(plot_UHD_data, county %in% c("Teton", "Canyon", "Twin Falls")), 
            inherit.aes = F,
            aes(x=ag_to_UHD+200, y=pvr_to_UHD-0.01, label=county), hjust="left") +
  geom_point(aes(x=ag_to_UHD, y=pvr_to_UHD, shape=barcode_R_lab, color = barcode_A_lab),
             data= plot_UHD_data, size=2) +
  labs(x="Ag Land Converted to Urban High Development (acres)",
       y="Avg PVR of Converted Land") +
  scale_color_viridis_d(name="Actor Type", option="turbo") +
  scale_shape_discrete(name = "Resource Type") +
  theme_bw()
# without R1-A1
ggplot(aes(x=ag_to_UHD, y=pvr_to_UHD, shape=barcode_R_lab, color = barcode_A_lab),
       data= subset(plot_UHD_data, barcode_A_lab!="A1")) +
  geom_hline(yintercept = mean_pvr_aoi, linetype="dashed") +
  geom_text(aes(x=2500, y=0.263, label="Avg Study Area PVR"), hjust="right", inherit.aes = F, data=data.frame()) +
  geom_segment(data=subset(plot_UHD_data, county %in% c("Teton", "Twin Falls")), 
               inherit.aes = F,
               aes(x=ag_to_UHD+200, xend=ag_to_UHD,
                   y=pvr_to_UHD-0.01, yend=pvr_to_UHD)) +
  geom_text(data=subset(plot_UHD_data, county %in% c("Teton", "Twin Falls")), 
            inherit.aes = F,
            aes(x=ag_to_UHD+200, y=pvr_to_UHD-0.01, label=county), hjust="left") +
  geom_point(aes(x=ag_to_UHD, y=pvr_to_UHD, shape=barcode_R_lab, color = barcode_A_lab),
             data= subset(plot_UHD_data, barcode_A_lab!="A1")) +
  labs(x="Ag Land Converted to Urban High Development (acres)",
       y="Avg PVR of Converted Land") +
  scale_color_viridis_d(name="Actor Type", option="turbo",
                        begin=0.2) +
  scale_shape_manual(name = "Resource Type",
                       values = c(17,15,3,7,8)) +
  theme_bw()

plot_LDR_data <- county_barcodes %>%
  rowwise() %>%
  mutate(ag_to_LDR = sum(c_across(c(crop_to_LDR, past_to_LDR, range_to_LDR)), na.rm=T),
         pvr_to_LDR = mean(c_across(c(crop_to_LDR_avg_pvr, past_to_LDR_avg_pvr, range_to_LDR_avg_pvr)), na.rm=T))

ggplot(aes(x=ag_to_LDR, y=pvr_to_LDR, shape=barcode_R_lab, color = barcode_A_lab),
       data= plot_LDR_data) +
  geom_hline(yintercept = mean_pvr_aoi, linetype="dashed") +
  geom_text(aes(x=12500, y=0.263, label="Avg Study Area PVR"), hjust="right", inherit.aes = F, data=data.frame()) +
  geom_segment(data=subset(plot_LDR_data, county %in% c("Teton", "Canyon", "Twin Falls")), 
               inherit.aes = F,
               aes(x=c(6300,850,2200), xend=ag_to_LDR,
                   y=c(0.345,0.2779,0.285), yend=pvr_to_LDR)) +
  geom_text(data=subset(plot_LDR_data, county %in% c("Teton", "Canyon", "Twin Falls")), 
            inherit.aes = F,
            aes(x=c(6300,700,2200), y=c(0.345,0.273,0.285), label=county), hjust="left") +
  geom_point(aes(x=ag_to_LDR, y=pvr_to_LDR, shape=barcode_R_lab, color = barcode_A_lab),
             data= plot_LDR_data, size=2) +
  labs(x="Ag Land Converted to Low Density Residential (acres)",
       y="Avg PVR of Converted Land") +
  scale_color_viridis_d(name="Actor Type", option="turbo") +
  scale_shape_discrete(name = "Resource Type") +
  theme_bw()
# without R1-A1
ggplot(aes(x=ag_to_LDR, y=pvr_to_LDR, shape=barcode_R_lab, color = barcode_A_lab),
       data= subset(plot_LDR_data, barcode_A_lab!="A1")) +
  geom_hline(yintercept = mean_pvr_aoi, linetype="dashed") +
  geom_text(aes(x=5500, y=0.283, label="Avg Study Area PVR"), hjust="right", inherit.aes = F, data=data.frame()) +
  geom_segment(data=subset(plot_LDR_data, county %in% c("Teton", "Twin Falls")), 
               inherit.aes = F,
               aes(x=c(850,2200), xend=ag_to_LDR,
                   y=c(0.2779,0.285), yend=pvr_to_LDR)) +
  geom_text(data=subset(plot_LDR_data, county %in% c("Teton", "Twin Falls")), 
            inherit.aes = F,
            aes(x=c(700,2200), y=c(0.273,0.285), label=county), hjust="left") +
  geom_point(aes(x=ag_to_LDR, y=pvr_to_LDR, shape=barcode_R_lab, color = barcode_A_lab),
             data= subset(plot_LDR_data, barcode_A_lab!="A1"), size=2) +
  labs(x="Ag Land Converted to Low Density Residential (acres)",
       y="Avg PVR of Converted Land") +
  scale_color_viridis_d(name="Actor Type", option="turbo",
                        begin=0.2) +
  scale_shape_manual(name = "Resource Type",
                     values = c(17,15,3,7,8)) +
  theme_bw()

##########################
# new AFT methods

nlcd01 <- FedData::get_nlcd(template = aoi, label="SRP",
                            year = 2001)
nlcd21 <- FedData::get_nlcd(template = aoi, label="SRP",
                            year = 2021)
imperv21 <- rast("/vsizip/R:/CarolynKoehn/Data_raw/Annual_NLCD_ImpDsc_2021_CU_C1V2.zip/Annual_NLCD_ImpDsc_2021_CU_C1V2.tif")
imperv21 <- crop(imperv21, st_transform(aoi, st_crs(imperv21)), mask=TRUE)
imperv21_prj <- project(imperv21, nlcd21, method="mode")

crops01 <- ifel(nlcd01==82, 1, NA)

grass01 <- ifel(nlcd01 %in% c(81, 71), 1, NA)

aft01 <- project(lulc_01, grass01, method = "mode")
grass01_aft <- ifel(grass01 == 1 & aft01 %in% c(2,3), 1, NA)

forest01 <- ifel(nlcd01 %in% c(41,42,43,90), 1, NA)

forest01_aft <- ifel(forest01 == 1 & aft01 == 5, 1, NA)

dev21 <- ifel(nlcd21 %in% c(21:24), 1, NA)

dev21_noroads <- ifel(imperv21_prj == 1, NA, dev21)

aglands <- ifel(crops01==1, 1, ifel(grass01_aft==1, 1, forest01_aft))
aglands_trim <- trim(aglands)
dev21_noroads_ext <- extend(trim(dev21_noroads), aglands_trim)
dev21_noroads_ext <- crop(dev21_noroads_ext, aglands_trim)
dev21_noroads_ext <- project(dev21_noroads_ext, aglands_trim, method="near")

converted <- ifel(aglands_trim==1 & dev21_noroads_ext == 1, 1, NA)

by_county <- extract(converted, st_transform(aoi, st_crs(converted)), fun=sum, na.rm=TRUE)

aglands_by_county <- extract(aglands, st_transform(aoi, st_crs(aglands)), fun=sum, na.rm=TRUE)

pvr_by_county <- project(pvr, converted, method="bilinear")
pvr_converted <- mask(pvr_by_county, converted)
pvr_bc2 <- extract(pvr_converted, st_transform(aoi, st_crs(pvr_by_county)), fun=mean, na.rm=TRUE)

aoi$convert_acres <- by_county$Class * 0.2223948433
aoi$ag_acres <- aglands_by_county$Class * 0.2223948433
aoi$convert_perc <- aoi$convert_acres / aoi$ag_acres * 100
aoi$convert_pvr <- pvr_bc2$productivity_versatility_resiliency_2016_conus_120m

county_barcodes <- left_join(county_barcodes, st_drop_geometry(aoi), by=c("county" = "NAME"))

# by acres
ggplot(aes(x=convert_acres, y=convert_pvr, shape=barcode_R_lab, color = barcode_A_lab),
       data= county_barcodes) +
  geom_hline(yintercept = mean_pvr_aoi, linetype="dashed") +
  geom_text(aes(x=15000, y=mean_pvr_aoi-0.01, label="Avg Study Area PVR"), hjust="right", inherit.aes = F, data=data.frame()) +
  geom_segment(data=subset(county_barcodes, county %in% c("Teton", "Canyon", "Twin Falls")), 
               inherit.aes = F,
               aes(x=convert_acres-50, xend=convert_acres,
                   y=convert_pvr-0.01, yend=convert_pvr)) +
  geom_text(data=subset(county_barcodes, county %in% c("Teton", "Canyon", "Twin Falls")), 
            inherit.aes = F,
            aes(x=convert_acres-50, y=convert_pvr-0.01, label=county), hjust="right") +
  geom_point(aes(x=convert_acres, y=convert_pvr, shape=barcode_R_lab, color = barcode_A_lab),
             data= county_barcodes, size=4) +
  labs(x="Ag Land Converted to Development (acres)",
       y="Avg PVR of Converted Land") +
  scale_color_viridis_d(name="Actor Type", option="turbo", direction=-1) +
  scale_shape_manual(name = "Resource Type", values = c(7,17,18,16,15,9)) +
  theme_bw()

# by percent converted
ggplot(aes(x=convert_perc, y=convert_pvr, shape=barcode_R_lab, color = barcode_A_lab),
       data= county_barcodes) +
  geom_hline(yintercept = mean_pvr_aoi, linetype="dashed") +
  geom_text(aes(x=15, y=mean_pvr_aoi-0.01, label="Avg Study Area PVR"), hjust="right", inherit.aes = F, data=data.frame()) +
  geom_segment(data=subset(county_barcodes, county %in% c("Teton", "Canyon", "Twin Falls")), 
               inherit.aes = F,
               aes(x=convert_perc+0.5, xend=convert_perc,
                   y=convert_pvr-0.01, yend=convert_pvr)) +
  geom_text(data=subset(county_barcodes, county %in% c("Teton", "Canyon", "Twin Falls")), 
            inherit.aes = F,
            aes(x=convert_perc+0.5, y=convert_pvr-0.01, label=county), hjust="left") +
  geom_point(aes(x=convert_perc, y=convert_pvr, shape=barcode_R_lab, color = barcode_A_lab),
             data= county_barcodes, size=4) +
  labs(x="Ag Land Converted to Development (%)",
       y="Avg PVR of Converted Land") +
  scale_color_viridis_d(name="Actor Type", option="turbo") +
  scale_shape_discrete(name = "Resource Type") +
  theme_bw()
