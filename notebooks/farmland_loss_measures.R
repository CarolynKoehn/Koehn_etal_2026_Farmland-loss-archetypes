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
                            levels = c("R1", "R5", "R4", "R2", "R3"),
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

# nlcd01 <- FedData::get_nlcd(template = aoi, label="SRP",
#                             year = 2001)
# nlcd21 <- FedData::get_nlcd(template = aoi, label="SRP",
#                             year = 2021)
nlcd01 <- rast("/vsizip/data/raw_data/Annual_NLCD_LndCov_2001_CU_C1V2.zip/Annual_NLCD_LndCov_2001_CU_C1V2.tif")
nlcd01 <- crop(nlcd01, st_transform(aoi, st_crs(nlcd01)), mask=TRUE)

nlcd25 <- rast("/vsizip/data/raw_data/Annual_NLCD_LndCov_2025_CU_C1V2.zip/Annual_NLCD_LndCov_2025_CU_C1V2.tif")
nlcd25 <- crop(nlcd25, st_transform(aoi, st_crs(nlcd25)), mask=TRUE)

# imperv21 <- rast("/vsizip/R:/CarolynKoehn/Data_raw/Annual_NLCD_ImpDsc_2021_CU_C1V2.zip/Annual_NLCD_ImpDsc_2021_CU_C1V2.tif")
# imperv21 <- crop(imperv21, st_transform(aoi, st_crs(imperv21)), mask=TRUE)
# imperv21_prj <- project(imperv21, nlcd21, method="mode")

imperv25 <- rast("/vsizip/data/raw_data/Annual_NLCD_ImpDsc_2025_CU_C1V2.zip/Annual_NLCD_ImpDsc_2025_CU_C1V2.tif")
imperv25 <- crop(imperv25, st_transform(aoi, st_crs(imperv25)), mask=TRUE)



crops01 <- ifel(nlcd01==82, 1, NA)

grass01 <- ifel(nlcd01 %in% c(81, 71), 1, NA)

aft01 <- project(lulc_01, grass01, method = "mode")
grass01_aft <- ifel(grass01 == 1 & aft01 %in% c(2,3,1102,1103), 1, NA)

forest01 <- ifel(nlcd01 %in% c(41,42,43,90), 1, NA)

forest01_aft <- ifel(forest01 == 1 & aft01 %in% c(5,1105), 1, NA)

dev25 <- ifel(nlcd25 %in% c(21:24), 1, NA)

dev25_noroads <- ifel(imperv25 == 1, NA, dev25)

aglands <- ifel(crops01==1, 1, ifel(grass01_aft==1, 1, forest01_aft))
# aglands_trim <- trim(aglands)
# dev21_noroads_ext <- extend(trim(dev21_noroads), aglands_trim)
# dev21_noroads_ext <- crop(dev21_noroads_ext, aglands_trim)
# dev21_noroads_ext <- project(dev21_noroads_ext, aglands_trim, method="near")

converted <- ifel(aglands==1 & dev25_noroads == 1, 1, NA)

by_county <- extract(converted, st_transform(aoi, st_crs(converted)), fun=sum, na.rm=TRUE)

aglands_by_county <- extract(aglands, st_transform(aoi, st_crs(aglands)), fun=sum, na.rm=TRUE)

pvr_prj <- project(pvr, converted, method="bilinear")
pvr_converted <- mask(pvr_prj, converted)
pvr_conv_by_county <- extract(pvr_converted, st_transform(aoi, st_crs(pvr_prj)), fun=mean, na.rm=TRUE)

aoi$convert_acres <- by_county$`NLCD Land Cover Class` * 0.2223948433
aoi$ag_acres <- aglands_by_county$`NLCD Land Cover Class` * 0.2223948433
aoi$convert_perc <- aoi$convert_acres / aoi$ag_acres * 100
aoi$convert_pvr <- pvr_conv_by_county$productivity_versatility_resiliency_2016_conus_120m

county_barcodes <- left_join(county_barcodes, st_drop_geometry(aoi), by=c("county" = "NAME"))

# by acres
ggplot(aes(x=convert_acres, y=convert_pvr, shape=barcode_R_lab, color = barcode_A_lab),
       data= county_barcodes) +
  geom_hline(yintercept = mean_pvr_aoi, linetype="dashed") +
  geom_text(aes(x=15000, y=mean_pvr_aoi-0.01, label="Avg Study Area PVR"), hjust="right", inherit.aes = F, data=data.frame()) +
  geom_segment(data=subset(county_barcodes, county %in% c("Teton", "Twin Falls")), 
               inherit.aes = F,
               aes(x=convert_acres+250, xend=convert_acres,
                   y=convert_pvr-0.01, yend=convert_pvr)) +
  geom_text(data=subset(county_barcodes, county %in% c("Teton", "Twin Falls")), 
            inherit.aes = F,
            aes(x=convert_acres+250, y=convert_pvr-0.01, label=county), hjust="left") +
  geom_segment(data=subset(county_barcodes, county == "Canyon"), 
               inherit.aes = F,
               aes(x=convert_acres-100, xend=convert_acres,
                   y=convert_pvr-0.01, yend=convert_pvr)) +
  geom_text(data=subset(county_barcodes, county == "Canyon"), 
            inherit.aes = F,
            aes(x=convert_acres-100, y=convert_pvr-0.01, label=county), hjust="right") +
  geom_point(aes(x=convert_acres, y=convert_pvr, shape=barcode_R_lab, color = barcode_A_lab),
             data= county_barcodes, size=4) +
  labs(x="Ag Land Converted to Development (acres)",
       y="Avg PVR of Converted Land") +
  scale_color_manual(name="Actor Type", 
                     values = c("black", viridis::turbo(n=6, direction=-1, begin=0.075, end=0.95)),
                     limits = c("A1", "A2", "A3", "A4➞A2", "A4", "A5➞A4" ,"A5" )) +
  #scale_color_viridis_d(name="Actor Type", option="turbo", direction=-1) +
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
  scale_color_manual(name="Actor Type", 
                     values = c("black", viridis::turbo(n=6, direction=-1, begin=0.075, end=0.95)),
                     limits = c("A1", "A2", "A3", "A4➞A2", "A4", "A5➞A4" ,"A5" )) +
  #scale_color_viridis_d(name="Actor Type", option="turbo") +
  scale_shape_manual(name = "Resource Type", values = c(7,17,18,16,15,9)) +
  theme_bw()

##########################
# use LCMAP as a different data source to verify?

lcmap90 <- rast("/vsizip/R:/CarolynKoehn/Data_raw/LCMAP/lcmap_df1bf009-933b-46cf-89e8-59b4f304ee21.zip/LCMAP_CU_1990_V13_LCPRI.tiff")
lcmap20 <- rast("/vsizip/R:/CarolynKoehn/Data_raw/LCMAP/lcmap_c632322b-8ba2-4342-bc85-46116c6b4dcb.zip/LCMAP_CU_2020_V13_LCPRI.tiff")

pvr <- rast("R:/SarahHalperin/Data/AFT_Data/productivity_versatility_resiliency_2016_conus_120m.tif")

aoi <- tigris::counties("ID", progress_bar = FALSE) %>%
  select(c("GEOID", "NAME")) %>%
  filter(NAME %in% c("Jefferson", "Gem", "Madison", "Clark", "Jerome", "Bonneville", "Teton", "Power", "Blaine", "Elmore", "Butte", "Adams", "Minidoka", "Bingham", "Ada", "Cassia", "Payette", "Fremont", "Bannock", "Twin Falls", "Owyhee", "Washington", "Gooding", "Lincoln", "Camas", "Canyon")) %>% st_transform(st_crs(lcmap20))

lulc90_srp <- crop(lcmap90, aoi, mask=TRUE)
lulc20_srp <- crop(lcmap20, aoi, mask=TRUE)
pvr_srp <- crop(project(pvr, lulc20_srp, method = "bilinear"), aoi, mask=TRUE)

crops90 <- ifel(lulc90_srp==2, 1, NA)
grass90 <- ifel(lulc90_srp==3, 1, NA)

# stopped here, can't check against AFT pastures/rangelands with these years

##############################

# farmland loss / developed gain by cluster

county_clusts_lossgain <- st_drop_geometry(clust_dat$Resource) %>%
  select(county, year, cluster, developed_change, Perc_change_land_in_farms, perc_farms_change) %>%
  left_join(st_drop_geometry(select(clust_dat$Actors, county, year, cluster)), 
            join_by(county, year), suffix = c(".R",".A")) %>%
  mutate(cluster.R = factor(cluster.R,
                            levels = c("R1", "R5", "R4", "R2", "R3"),
                            labels = c("R1", "R2", "R3", "R4", "R5")),
         cluster.A = factor(cluster.A,
                            levels = c("A1", "A2", "A5", "A3", "A4"),
                            labels = c("A1", "A2", "A3", "A4", "A5")))

ggplot(county_clusts_lossgain, aes(x=developed_change, y=perc_farms_change,
                                   color = cluster.A, shape = cluster.R)) +
  geom_point(size=3)
ggplot(county_clusts_lossgain, aes(x=developed_change, y=Perc_change_land_in_farms,
                                   color = cluster.A, shape = cluster.R)) +
  geom_point(size=3)
# not easy to interpret

##############################

# AFT methods with NLCD - but get a slope instead of two points in time

nlcd01 <- rast("/vsizip/data/raw_data/Annual_NLCD_LndCov_2001_CU_C1V2.zip/Annual_NLCD_LndCov_2001_CU_C1V2.tif")
nlcd01 <- crop(nlcd01, st_transform(aoi, st_crs(nlcd01)), mask=TRUE)

nlcd_05_25 <- list(rast(paste0("/vsizip/data/raw_data/NLCD_LndCov_1985-2020.zip/",
                            unzip("data/raw_data/NLCD_LndCov_1985-2020.zip", list=TRUE)$Name[str_detect(unzip("data/raw_data/NLCD_LndCov_1985-2020.zip", list=TRUE)$Name, ".tiff$")][c(6:8,1)])), 
                     rast("/vsizip/data/raw_data/Annual_NLCD_LndCov_2025_CU_C1V2.zip/Annual_NLCD_LndCov_2025_CU_C1V2.tif"))
nlcd_05_25_crop <- rast(lapply(nlcd_05_25, function(r) 
  crop(r, st_transform(aoi, st_crs(nlcd_05_25[[1]])), mask=TRUE)))
rm(nlcd_05_25)

imperv_05_25 <- rast(paste0("/vsizip/data/raw_data/Annual_NLCD_ImpDsc_", seq(2005,2025,5), "_CU_C1V2.zip/Annual_NLCD_ImpDsc_", seq(2005,2025,5), "_CU_C1V2.tif"))
imperv_05_25_crop <- crop(imperv_05_25, st_transform(aoi, st_crs(imperv_05_25)), mask=TRUE)
rm(imperv_05_25)

crops01 <- ifel(nlcd01==82, 1, NA)

grass01 <- ifel(nlcd01 %in% c(81, 71), 1, NA)

aft01 <- project(lulc_01, grass01, method = "mode")
grass01_aft <- ifel(grass01 == 1 & aft01 %in% c(2,3,1102,1103), 1, NA)

forest01 <- ifel(nlcd01 %in% c(41,42,43,90), 1, NA)

forest01_aft <- ifel(forest01 == 1 & aft01 %in% c(5,1105), 1, NA)

dev_05_25 <- ifel(nlcd_05_25_crop %in% c(21:24), 1, NA)

dev_05_25_noroads <- dev_05_25
for(i in 1:nlyr(dev_05_25_noroads)) {
  dev_05_25_noroads[[i]] <- ifel(imperv_05_25_crop[[i]] == 1, NA, dev_05_25_noroads[[i]])
}

aglands <- ifel(crops01==1, 1, ifel(grass01_aft==1, 1, forest01_aft))

converted <- dev_05_25_noroads
for(i in 1:nlyr(converted)) {
  converted[[i]] <- ifel(aglands==1 & dev_05_25_noroads[[i]] == 1, 1, NA)
}
names(converted) <- paste0("converted_", seq(2005,2025,5))

by_county <- extract(converted, st_transform(aoi, st_crs(converted)), fun=sum, na.rm=TRUE)

aglands_by_county <- extract(aglands, st_transform(aoi, st_crs(aglands)), fun=sum, na.rm=TRUE)

pvr_prj <- project(pvr, aglands, method="bilinear")
pvr_converted_25 <- mask(pvr_prj, converted$converted_2025)
pvr_conv_by_county <- extract(pvr_converted_25, st_transform(aoi, st_crs(pvr_prj)), fun=mean, na.rm=TRUE)

by_county_rates <- cbind(by_county, aoi) %>%
  cbind(., pvr_conv_by_county$productivity_versatility_resiliency_2016_conus_120m) %>%
  rename(converted_pvr = `pvr_conv_by_county$productivity_versatility_resiliency_2016_conus_120m`) %>%
  cbind(., aglands_by_county$`NLCD Land Cover Class`) %>%
  rename("aglands01" = "aglands_by_county$`NLCD Land Cover Class`") %>%
  st_drop_geometry() %>%
  select(-c(ID, GEOID, geometry)) %>%
  mutate(converted_2001=0) %>%
  pivot_longer(!c(NAME, converted_pvr, aglands01), 
               names_to = "year",
               values_to = "converted_px") %>%
  mutate(year = as.integer(str_extract(year, "[0-9]+")),
         converted_px = ifelse(is.nan(converted_px), 0, converted_px),
         converted_acres = converted_px * 0.2223948433,
         aglands01_acres = aglands01 * 0.2223948433)

county_barcodes[,c("rate", "CIlow", "CIhigh")] <- NA
for(i in 1:nrow(county_barcodes)) {
  county_barcodes$rate[i] <- 
    summary(lm(converted_acres~year, 
               data = subset(by_county_rates, 
                             NAME==county_barcodes$county[i])))$coefficients[2,1]
  
  county_barcodes$CIlow[i] <- 
    confint(lm(converted_acres~year, 
               data = subset(by_county_rates, 
                             NAME==county_barcodes$county[i])))[2,1]
  
  county_barcodes$CIhigh[i] <- 
    confint(lm(converted_acres~year, 
               data = subset(by_county_rates, 
                             NAME==county_barcodes$county[i])))[2,2]
}

county_barcodes <- county_barcodes %>%
  left_join(unique(by_county_rates[,c("NAME", "converted_pvr", "aglands01_acres")]), 
            join_by(county == NAME))

county_barcodes_percCIs <- county_barcodes %>%
  mutate(perc_rate = rate / aglands01_acres * 100,
         perc_rate_CIlow = CIlow / aglands01_acres * 100,
         perc_rate_CIhigh = CIhigh / aglands01_acres * 100)


# acres per year
ggplot(aes(x=rate, y=converted_pvr, shape=barcode_R_lab, color = barcode_A_lab),
       data= county_barcodes) +
  geom_hline(yintercept = mean_pvr_aoi, linetype="dashed") +
  geom_text(aes(x=15, y=mean_pvr_aoi-0.01, label="Avg Study Area PVR"), hjust="right", inherit.aes = F, data=data.frame()) +
  geom_segment(data=subset(county_barcodes, county %in% c("Teton", "Canyon", "Twin Falls")), 
               inherit.aes = F,
               aes(x=rate+0.5, xend=rate,
                   y=converted_pvr-0.01, yend=converted_pvr)) +
  geom_point(aes(x=rate, y=converted_pvr, shape=barcode_R_lab, color = barcode_A_lab),
             data= county_barcodes, size=4) +
  geom_text(data=subset(county_barcodes, county %in% c("Teton", "Canyon", "Twin Falls")), 
            inherit.aes = F,
            aes(x=rate+0.5, y=converted_pvr-0.01, label=county), hjust="left") +
  labs(x="Ag Land Converted to Development (acres/yr)",
       y="Avg PVR of Converted Land") +
  scale_color_manual(name="Actor Type", 
                     values = c("black", viridis::turbo(n=6, direction=-1, begin=0.075, end=0.95)),
                     limits = c("A1", "A2", "A3", "A4➞A2", "A4", "A5➞A4" ,"A5" )) +
  scale_shape_manual(name = "Resource Type", values = c(7,17,18,16,15,9),
                     limits = c("R1", "R2", "R3➞R2", "R3", "R4" ,"R5" )) +
  theme_bw()

# percent change per year
# conclusion: the results are not different enough to justify the extra compelixty -- there don't seem to be anomalous 2025 measurements skewing where the average rate wouldn't go
ggplot(aes(x=perc_rate, y=converted_pvr, shape=barcode_R_lab, color = barcode_A_lab),
       data= county_barcodes_percCIs) +
  geom_hline(yintercept = mean_pvr_aoi, linetype="dashed") +
  geom_text(aes(x=1, y=mean_pvr_aoi-0.01, label="Avg Study Area PVR"), hjust="right", inherit.aes = F, data=data.frame()) +
  geom_segment(data=subset(county_barcodes_percCIs, county %in% c("Teton", "Canyon", "Twin Falls")), 
               inherit.aes = F,
               aes(x=perc_rate+0.5, xend=perc_rate,
                   y=converted_pvr-0.01, yend=converted_pvr)) +
  geom_point(aes(x=perc_rate, y=converted_pvr, shape=barcode_R_lab, color = barcode_A_lab),
             data= county_barcodes_percCIs, size=4) +
  geom_text(data=subset(county_barcodes_percCIs, county %in% c("Teton", "Canyon", "Twin Falls")), 
            inherit.aes = F,
            aes(x=perc_rate+0.5, y=converted_pvr-0.01, label=county), hjust="left") +
  labs(x="Annual Percent Conversion of Ag Land to Development",
       y="Avg PVR of Converted Land") +
  scale_color_manual(name="Actor Type", 
                     values = c("black", viridis::turbo(n=6, direction=-1, begin=0.075, end=0.95)),
                     limits = c("A1", "A2", "A3", "A4➞A2", "A4", "A5➞A4" ,"A5" )) +
  scale_shape_manual(name = "Resource Type", values = c(7,17,18,16,15,9),
                     limits = c("R1", "R2", "R3➞R2", "R3", "R4" ,"R5" )) +
  theme_bw()
