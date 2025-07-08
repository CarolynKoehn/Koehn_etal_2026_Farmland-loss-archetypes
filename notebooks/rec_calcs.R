library(sf)
library(terra)
library(tidyverse)
library(exactextractr)

nps <- st_read("R:/Common/CarolynKoehn/Data_raw/Parks and Recreation/NPS_-_Land_Resources_Division_Boundary_and_Tract_Data_Service/NPS_-_Land_Resources_Division_Boundary_and_Tract_Data_Service.shp")

nps_id <- subset(nps, STATE == "ID")

idaho <- tigris::states() %>%
  filter(STUSPS == "ID") %>%
  st_transform("epsg:5070")

pad_id <- st_read("R:/CarolynKoehn/Data_raw/PAD-US 3.0/PADUS3_0_State_ID_SHP/PADUS3_0Combined_StateID.shp") %>%
  st_transform("epsg:5070") %>%
  st_make_valid()

pad_rec <- pad_id %>%
  filter(d_Des_Tp %in% c("Recreation Management Area", "National Monument", "National Park", "State Recreation Area", "State Park", "National Recreation Area"))

aoi <- st_read("R:/Common/CarolynKoehn/Data_clipped/00_SnakeRiverPlainCounties_StudyArea/SRP_counties_26.shp") %>%
  st_transform("epsg:5070")

pad_rec_clip <- pad_rec %>%
  st_transform("epsg:5070") %>%
  st_intersection(., idaho)

#template
ext_template <- rast(vect(idaho))
res(ext_template) <- 10

pad_rec_years <- vector("list", length=4)
names(pad_rec_years) <- seq(1990, 2020, 10)

pad_rec_years <- mapply(function(list, year) {
  pad_rec_clip %>% filter(Date_Est <= year | is.na(Date_Est))
}, list = pad_rec_years, year = seq(1990, 2020, 10), SIMPLIFY = FALSE)

pad_rec_r <- lapply(pad_rec_years, function(df) {
  rasterize(df,
            ext_template, field = 1)
})

pad_rec_r <- rast(pad_rec_r)
  
sum_cover <- function(x){
  all_props <- list(x %>%
         group_by(value) %>%
         summarize(total_area = sum(coverage_area)) %>%
         mutate(proportion = total_area/sum(total_area)))
  
  return(subset(all_props[[1]], value == 1)$proportion * 100)
  
}

#extract the area of each raster cell covered by the plot and summarize
x <- exact_extract(pad_rec_r, aoi, 
                   coverage_area = TRUE, summarize_df = TRUE, fun = sum_cover,
                   stack_apply = TRUE, include_cols = "NAME")


#add plot names to the elements of the output list
rownames(x) <- aoi$NAME


######################
library(tidyverse)
library(tidycensus)
library(rnassqs)
library(terra)
library(sf)
library(FedData)
library(exactextractr)
library(censusapi)

# set API keys
readRenviron("../.Renviron")
nassqs_auth(key = Sys.getenv('NASSQS_TOKEN'))
census_api_key(Sys.getenv('CENSUS_API_KEY'))

# county data
county_shapes <- tigris::counties("ID") %>%
  select(c("GEOID", "NAME")) %>%
  mutate(tot_area = st_area(.))

# names and GEOIDs
GEOID_names <- st_drop_geometry(county_shapes$NAME)
names(GEOID_names) <- st_drop_geometry(county_shapes$GEOID)


# travel time to
# doesn't work right now
id_roads <- tigris::roads(state=16,
                           county=substr(county_shapes$GEOID, 3, 5),
                           year=2020)
wy_roads <- tigris::roads(state = "WY",
                          county = "Park", "Teton",
                          year = 2020)
mo_roads <- tigris::roads(state = "MT",
                          county = "Beaverhead", "Madison", "Gallatin",
                          year = 2020)

all_roads_car <- rbind(id_roads, wy_roads, mo_roads) %>%
  filter(MTFCC %in% c("S1100", "S1200", "S1400", "S1630")) %>%
  mutate(road_cat = factor(MTFCC, 
                        levels = c("S1100", "S1200", "S1400"),
                        labels = c("Interstate", "Highway", "Road")),
         speed = ifelse(road_cat == "Interstate", 80,
                        ifelse(road_cat == "Highway", 60, 25)))

roads.r <- rasterize(all_roads_car, ext_template, field="speed", fun="max")
