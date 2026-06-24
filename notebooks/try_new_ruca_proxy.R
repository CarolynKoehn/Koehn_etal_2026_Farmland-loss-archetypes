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

############

# by urban populations

# retrieve all US counties
us_counties <- tigris::counties()
# create 150 mile buffer around Idaho
id_buff <- st_union(county_shapes) %>% st_buffer(units::set_units(150, "mi"))
# use only counties falling within buffer
near_counties <- st_filter(us_counties, id_buff, .predicate = st_covered_by)

###### get block data from urban land code chunk

# for each year of analysis:
years <- seq(1990, 2020, 10)
# need the block data list from urban areas code chunk, but only the first version (do not need distance from urban labels)
counties_urb <- mapply(function(blocks, year) {
  # read in urban areas data saved from the urban land code chunk - RS3
  UAs <- st_read(paste0("data/urb_areas_thresholded/urb_areas_thresholded_", year, ".shp"))
  
  # retain blocks in urban areas
  UA_blocks <- st_join(blocks, UAs, join = st_covered_by)
  
  # get population of each urban area to categorize as metropolitan, micropolitan, small town
  UA_pops <- UA_blocks %>%
    st_drop_geometry() %>%
    group_by(FID) %>%
    summarise(UA_pop = sum(tot_pop)) %>%
    na.omit() %>%
    right_join(UAs) %>%
    st_as_sf() %>%
    mutate(urb_label = ifelse(UA_pop >= 50000 & !is.na(UA_pop), "Metropolitan", NA)) %>%
    mutate(urb_label = ifelse(UA_pop >= 10000 & UA_pop < 50000 & !is.na(UA_pop), "Micropolitan", urb_label)) %>%
    # smallest UAs classified as "small town"
    mutate(urb_label = ifelse(is.na(urb_label), "Small town", urb_label))
  
  # determine which urban areas are in each county
  county_UAs <- st_intersects(us_counties, st_transform(UA_pops, st_crs(us_counties)))
  # find the largest urban area in each county
  county_UAs_max_pop <- unlist(lapply(county_UAs, function(x) {
    if(length(x) == 0) {
      return_val <- NA
    } else {
      pops <- st_drop_geometry(UA_pops[x, "UA_pop"])
      return_val <- max(pops)
    }
    return(return_val)
  }))
  
  # label counties by their largest urban area population
  us_counties$max_urb_pop <- county_UAs_max_pop
  
  # add RUCA-style labels based on population
  county_UAs_labels <- us_counties %>%
    mutate(urb_label = ifelse(!(GEOID %in% near_counties$GEOID), ">150 mi commute", NA)) %>%
    mutate(urb_label = ifelse(max_urb_pop >= 50000 & !is.na(max_urb_pop), "Metropolitan", urb_label)) %>%
    mutate(urb_label = ifelse(max_urb_pop >= 10000 & max_urb_pop < 50000 & !is.na(max_urb_pop), "Micropolitan", urb_label)) %>%
    # all small urban areas = "small town"
    mutate(urb_label = ifelse(!is.na(max_urb_pop) & is.na(urb_label), "Small town", urb_label)) %>%
    mutate(urb_label = ifelse(is.na(urb_label), "Rural", urb_label))
  
  return(county_UAs_labels)
}, blocks = block_data_list, year = years, SIMPLIFY = FALSE)

# calculate block level populations in urban categories for each county
id_counties_urb_pops <- mapply(function(blocks, year) {
  # read in urban areas data saved from the urban land code chunk - RS3
  UAs <- st_read(paste0("data/urb_areas_thresholded/urb_areas_thresholded_", year, ".shp"))
  
  # retain blocks in urban areas
  UA_blocks <- st_join(blocks, UAs, join = st_covered_by)
  
  # get population of each urban area to categorize as metropolitan, micropolitan, small town
  UA_pops <- UA_blocks %>%
    st_drop_geometry() %>%
    group_by(FID) %>%
    summarise(UA_pop = sum(tot_pop)) %>%
    na.omit() %>%
    right_join(UAs) %>%
    st_as_sf() %>%
    mutate(urb_label = ifelse(UA_pop >= 50000 & !is.na(UA_pop), "Metropolitan", NA)) %>%
    mutate(urb_label = ifelse(UA_pop >= 10000 & UA_pop < 50000 & !is.na(UA_pop), "Micropolitan", urb_label)) %>%
    # smallest UAs classified as "small town"
    mutate(urb_label = ifelse(is.na(urb_label), "Small town", urb_label))
  
  # add urban labels to blocks
  blocks_urblabs <- blocks %>%
    left_join(st_drop_geometry(UA_blocks[,c("GISJOIN", "FID")]), by="GISJOIN") %>%
    left_join(st_drop_geometry(UA_pops[,c("FID", "urb_label")]), by="FID") %>%
    mutate(urb_label = ifelse(is.na(urb_label), "Rural", urb_label))
  
  # total populations by urban label and county
  urb_pops_cty <- blocks_urblabs %>%
    st_drop_geometry(.) %>%
    group_by(FIPS, urb_label) %>%
    summarise(urb_pop = sum(tot_pop, na.rm=TRUE)) %>%
    pivot_wider(names_from = urb_label,
                values_from = urb_pop) %>%
    ungroup() %>%
    replace(is.na(.), 0) %>%
    rowwise(FIPS) %>%
    mutate(tot_pop = sum(c_across(everything()))) %>%
    mutate(across(everything(), ~.x/tot_pop * 100)) %>%
    select(-tot_pop) %>%
    filter(str_starts(FIPS, "16")) %>%
    mutate(year = year)
  
  return(urb_pops_cty)
}, blocks = block_data_list, year = years, SIMPLIFY = FALSE)

# get commuting flows
commutes_1990 <- read.table("https://www2.census.gov/programs-surveys/commuting/datasets/1990/worker-flow/idresco.txt", sep = "\t", strip.white = TRUE) %>%
  # format data with string splitting based on metadata
  mutate(allinfo = str_squish(V1),
         ResFIPS = str_replace_all(str_sub(allinfo, 1, 6), " ", ""),
         WorkFIPS = str_replace_all(str_sub(allinfo, 14, 19), " ", ""),
         Workers = as.integer(str_extract(str_extract(allinfo, " \\d+ [A-Z]"), "\\d+"))) %>%
  # join 1990 county urbanicity data for resident counties
  left_join(st_drop_geometry(counties_urb[[1]])[,c("GEOID", "NAME", "urb_label", "max_urb_pop")], by = c("ResFIPS" = "GEOID")) %>%
  rename("ResUrbStatus" = "urb_label",
         "ResUrbPop" = "max_urb_pop",
         "county" = "NAME") %>%
  # join 1990 county urbanicity data for work counties
  left_join(st_drop_geometry(counties_urb[[1]])[,c("GEOID", "urb_label", "max_urb_pop")], by = c("WorkFIPS" = "GEOID")) %>%
  rename("WorkUrbStatus" = "urb_label",
         "WorkUrbPop" = "max_urb_pop") %>%
  rowwise() %>%
  # label commutes based on resident and work county
  mutate(WorkUrbStatus = ifelse(is.na(WorkUrbStatus), ">150 mi commute", WorkUrbStatus),
         commute = ifelse(ResFIPS == WorkFIPS, 
                          "Stay",
                          ifelse(ResUrbStatus == "Metropolitan" & WorkUrbStatus == "Metropolitan",
                                 ifelse(WorkUrbPop > ResUrbPop, 
                                        "BiggerMetropolitan",
                                        "SmallerMetropolitan"),
                                 WorkUrbStatus)
         ),
         year = 1990) %>%
  group_by(county, year, commute, ResFIPS) %>%
  summarise(workers = sum(Workers))

commutes_2000 <- readxl::read_xls("data/raw_data/2kresco_id.xls",
                                  # leave out empty rows
                                  range = "A5:K1388") %>%
  rowwise() %>%
  # format data with string splitting based on metadata
  mutate(ResFIPS = paste0(`Res State`, `Res County`),
         WorkFIPS = paste0(str_sub(`Work State`, 2,3), `Work County`)) %>%
  # join 2000 county urbanicity data for resident counties
  left_join(st_drop_geometry(counties_urb[[2]])[,c("GEOID", "NAME", "urb_label", "max_urb_pop")], by = c("ResFIPS" = "GEOID")) %>%
  rename("ResUrbStatus" = "urb_label",
         "ResUrbPop" = "max_urb_pop",
         "county" = "NAME") %>%
  # join 2000 county urbanicity data for work counties
  left_join(st_drop_geometry(counties_urb[[2]])[,c("GEOID", "urb_label", "max_urb_pop")], by = c("WorkFIPS" = "GEOID")) %>%
  rename("WorkUrbStatus" = "urb_label",
         "WorkUrbPop" = "max_urb_pop") %>%
  # label commutes based on resident and work county
  mutate(WorkUrbStatus = ifelse(is.na(WorkUrbStatus), ">150 mi commute", WorkUrbStatus),
         commute = ifelse(ResFIPS == WorkFIPS, 
                          "Stay",
                          ifelse(ResUrbStatus == "Metropolitan" & WorkUrbStatus == "Metropolitan",
                                 ifelse(WorkUrbPop > ResUrbPop, 
                                        "BiggerMetropolitan",
                                        "SmallerMetropolitan"),
                                 WorkUrbStatus)
         ),
         year = 2000) %>%
  group_by(county, year, commute, ResFIPS) %>%
  summarise(workers = sum(Count))

commutes_2010 <- readxl::read_xlsx("data/raw_data/table1 (2010).xlsx",
                                   # leave out empty rows, undesired column names
                                   range = "A6:J26342",
                                   # manually name columns
                                   col_names = c("Res_State_FIPS", "Res_County_FIPS", "Work_State_FIPS", "Work_County_FIPS", "Workers", "MOE", "Res_State", "Res_County", "Work_State", "Work_County")) %>%
  # filter to Idaho
  filter(Res_State == "Idaho") %>%
  rowwise() %>%
  mutate(ResFIPS = paste0(Res_State_FIPS, Res_County_FIPS),
         WorkFIPS = paste0(str_sub(Work_State_FIPS, 2,3), Work_County_FIPS)) %>%
  # join 2010 county urbanicity data for resident counties
  left_join(st_drop_geometry(counties_urb[[3]])[,c("GEOID", "NAME", "urb_label", "max_urb_pop")], by = c("ResFIPS" = "GEOID")) %>%
  rename("ResUrbStatus" = "urb_label",
         "ResUrbPop" = "max_urb_pop",
         "county" = "NAME") %>%
  # join 2010 county urbanicity data for work counties
  left_join(st_drop_geometry(counties_urb[[3]])[,c("GEOID", "urb_label", "max_urb_pop")], by = c("WorkFIPS" = "GEOID")) %>%
  rename("WorkUrbStatus" = "urb_label",
         "WorkUrbPop" = "max_urb_pop") %>%
  # label commutes based on resident and work county
  mutate(WorkUrbStatus = ifelse(is.na(WorkUrbStatus), ">150 mi commute", WorkUrbStatus),
         commute = ifelse(ResFIPS == WorkFIPS, 
                          "Stay",
                          ifelse(ResUrbStatus == "Metropolitan" & WorkUrbStatus == "Metropolitan",
                                 ifelse(WorkUrbPop > ResUrbPop, 
                                        "BiggerMetropolitan",
                                        "SmallerMetropolitan"),
                                 WorkUrbStatus)
         ),
         year = 2010) %>%
  group_by(county, year, commute, ResFIPS) %>%
  summarise(workers = sum(Workers))

commutes_2020 <- readxl::read_xlsx("data/raw_data/table1 (2020).xlsx",
                                   # leave out empty rows, undesired column names
                                   range = "A9:J23788",
                                   # manually rename columns
                                   col_names = c("Res_State_FIPS", "Res_County_FIPS", "Res_State", "Res_County", "Work_State_FIPS", "Work_County_FIPS", "Work_State", "Work_County", "Workers", "MOE" )) %>%
  # filter to Idaho
  filter(Res_State == "Idaho") %>%
  rowwise() %>%
  mutate(ResFIPS = paste0(Res_State_FIPS, Res_County_FIPS),
         WorkFIPS = paste0(str_sub(Work_State_FIPS, 2,3), Work_County_FIPS)) %>%
  # join 2020 county urbanicity data for resident counties
  left_join(st_drop_geometry(counties_urb[[4]])[,c("GEOID", "NAME", "urb_label", "max_urb_pop")], by = c("ResFIPS" = "GEOID")) %>%
  rename("ResUrbStatus" = "urb_label",
         "ResUrbPop" = "max_urb_pop",
         "county" = "NAME") %>%
  # join 2020 county urbanicity data for work counties
  left_join(st_drop_geometry(counties_urb[[4]])[,c("GEOID", "urb_label", "max_urb_pop")], by = c("WorkFIPS" = "GEOID")) %>%
  rename("WorkUrbStatus" = "urb_label",
         "WorkUrbPop" = "max_urb_pop") %>%
  # label commutes based on resident and work county
  mutate(WorkUrbStatus = ifelse(is.na(WorkUrbStatus), ">150 mi commute", WorkUrbStatus),
         commute = ifelse(ResFIPS == WorkFIPS, 
                          "Stay",
                          ifelse(ResUrbStatus == "Metropolitan" & WorkUrbStatus == "Metropolitan",
                                 ifelse(WorkUrbPop > ResUrbPop, 
                                        "BiggerMetropolitan",
                                        "SmallerMetropolitan"),
                                 WorkUrbStatus)
         ),
         year = 2020) %>%
  group_by(county, year, commute, ResFIPS) %>%
  summarise(workers = sum(Workers))

# combine data from all years
commutes_all <- do.call(rbind, list(commutes_1990, commutes_2000, commutes_2010, commutes_2020)) %>%
  group_by(county, year) %>%
  # find total workers for each county-year
  mutate(tot_workers = sum(workers)) %>%
  rowwise() %>%
  # percent workers with each type of commute in each county-year
  mutate(perc_workers = workers / tot_workers * 100) %>%
  # pivot so that each commute type is a column, fill NAs with 0
  pivot_wider(id_cols = c("county", "year", "ResFIPS"),
              names_from = "commute",
              names_prefix = "commute_",
              values_from = "perc_workers",
              values_fill = 0)

urbs_plus_commutes <- do.call("rbind", id_counties_urb_pops) %>%
  right_join(commutes_all, by = c("FIPS" = "ResFIPS", "year")) %>%
  relocate(county, year, .before=FIPS) %>%
  ungroup() %>%
  replace(is.na(.), 0) %>%
  select(-FIPS)

ruca_proxy_pca <- prcomp(urbs_plus_commutes[,-c(1:2)])
summary(ruca_proxy_pca)
