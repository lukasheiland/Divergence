
# Notes -------------------------------------------------------------------

## Actual subsetting protocol
# 0. Observations were dropped with NA environmental variables.
# 1. Within species, only survey iterations per plot were kept, where both life stages of the species had been assessed in the survey, so that the sampled environmental space was constant across life stages and surveys.
# 2. Subsequently, within country and species, data were subset to include the first and the last survey (to coherently have two survey iterations per country to make estimates comparable).
# 3. Finally, species with less than 20 observations per country, life stage and survey iteration have been excluded.
    # (3. is done in the Fit* script.)

# The Library -------------------------------------------------------------
library(rprojroot)

library(tidyverse)
library(magrittr)
library(scales)
library(sf)

library(lubridate)

library(missRanger)
library(parallel)
n_cores <-  getOption("mc.cores", 3L)

# Set working directory ---------------------------------------------------
# Set the working directory to the git root, which is presumed to be the parent directory of the source file location
setwd(rprojroot::find_root(rprojroot::criteria$is_git_root))

# Source ------------------------------------------------------------------
# this references the forest inventory data directory, which is  another project
inventoriespath <- '../Inventories/' 

## source the data preprocessing script:
# source(paste0(inventoriespath, 'Master.R'))


# Seed ------------------------------------------------------------------
# Set seed, here for random forest imputation reproducibility.
set.seed(1)


# !!! -----------------------------------------------------------------------
## Will be used for species selection and plot selection further down.
obsthreshold <- 20 # minimal n_obs within species/inventory/obsid/ontogeny


# Predictors --------------------------------------------------------------
Env <- readRDS(paste(inventoriespath, 'Combined data', 'Env.rds', sep = '/'))
Env$lat <-  st_coordinates(Env)[,2]

predictor_select <- c("cwbYear_aclim",
                      # "cwbGdd0_aclim",
                      "gdd0_aclim" #,
                      # "tPeriodic2010_mh",
                      # "precPeriodic2010_mh",
                      # "clay_esdact",
                      # "sand_esdact",
                      # "octop_esdacoc",
                      # "cn_esdacc",
                      # "p_esdacc",
                      # "phCaCl_esdacc",
                      # "tMinColdMonth_wc2",
                      # "wwpi_cop"
                      # "lat"
                      )

idname <- c("envjoinid", "plotid", "plotobsid", "obsid", "clusterobsid", "clusterid", "time")

Env %<>% dplyr::select(all_of(c(idname, predictor_select)))

## plotid is not necessarily the id for unique environments (plotid == envjoinid as it is DE_BWI).
## As is, the predictors are available on the plotobsid level in SE_NFI. But here Env is only needed for plotids.
## The coordinates do not change per observation of the same plot (extracted variables are dependent on coordinates), but the local variables aspect_loc and slope_loc do.

Env %<>% dplyr::filter(!duplicated(plotid)) #!!! It's ok in this case. But only if environment is joined by plotid (see below).

## Save before NA handling.
saveRDS(Env, 'Data/Env.rds')

## Count NAs
Env %>% is.na() %>% colSums() # plotobsid/time are purposefully NA in DE_BWI as the environment does only change with plotid, i.e. is constant over times (time).

## Explore NA locations
require(leaflet)
basemap <- leaflet(Env) %>% fitBounds(3, 62, 10, 40) %>% addTiles(urlTemplate = 'https://{s}.tile.opentopomap.org/{z}/{x}/{y}.png')
nacoord <- st_coordinates(Env[is.na(Env$gdd0_aclim), "gdd0_aclim"])
addCircles(basemap, lng = nacoord[,1], lat = nacoord[,2], color = 'red') # NAs in most datasets are at the coasts, on land close to water.
## _esdacc NA are mostly close to water.
## _wc2 NA for coastline, etc.
## _cgiar NA for coastline
## _esdact NA close to coastlines and to Swiss border
## _mh NA for coastline

## IMPUTATION.
## Problem: Water-logged places often are NA in rasters.
## Solution: those should be easily predictable by features and coordinates.
imputeSpatial <-  function(Sf, varname_lhs, varname_rhs, ...){
  imputationformula <- paste(paste(varname_lhs, collapse = '+'), '~ X + Y +', paste(varname_rhs, collapse = '+')) %>%
    as.formula()
  Sf %<>% bind_cols(as.data.frame(st_coordinates(.)))
  geom <<- st_geometry(Sf)
  Sf %<>% st_drop_geometry() %>%
    missRanger(imputationformula, ...) %>%
    dplyr::select(-X, -Y) %>%
    st_set_geometry(geom)
  return(Sf)
}

varname_lhs <- predictor_select # The variables to be imputed. All variables have waterfront NAs.
varname_rhs <- predictor_select # The features to be predicted from. X, Y coordinates will be added by imputeSpatial().

Env_imputed <- imputeSpatial(Env, varname_lhs,
                             varname_rhs,
                             returnOOB = T,
                             seed = 1,
                             num.trees = 25)

## So many NAs
Env %>% is.na() %>% colSums() # plotobsid/time are purposefully NA in DE_BWI (Env changes only on plot level, not in time.)
## of
varname_lhs
## have been imputed with
varname_rhs
##. Has taken 4 iterations with num.trees = 25.
## With final average out of bag prediction error (oob):
# dput(attr(Env_imputed, "oob"))
# c(cwbYear_aclim = 0.000703498097850094, cwbGdd0_aclim = 0.000839791244127095, 
#   gdd0_aclim = 0.00070294745930576, tMinColdMonth_wc2 = 0.0006023741105646, 
#   wwpi_cop = 0.012081488887649, clay_esdact = 0.00190083851873737, 
#   sand_esdact = 0.00160083954955943, p_esdacc = 0.00581990457037878, 
#   cn_esdacc = 0.0019057081428244, phCaCl_esdacc = 0.00253252059989526)

## Save version with imputation.
saveRDS(Env_imputed, 'Data/Env_imputed.rds')

## IMPUTED DATA IS NOT UTILIZED!
## !!!
# Env <- Env_imputed
rm(Env_imputed)

# Load data ---------------------------------------------------------------
Small_pres <- readRDS(paste0(inventoriespath, 'Combined data/Small_pres.rds')) # this is long and tidy and zero-completed.
Big_pres <- readRDS(paste0(inventoriespath, 'Combined data/Big_pres.rds'))


# Select species ----------------------------------------------------------
## Species are filtered here, to improve wrangling speed,
## but the same criteria will be applied again later, after NA dropping etc.

Taxa <- read.csv(paste0(inventoriespath, 'Taxa/Taxa.csv'))

tax <- unique(as.character(Taxa$tax))
species <- tax[grepl("^[^.]*\\.[^.]*$", tax)] # select those taxa that are species

Species <- bind_rows(Small_pres, Big_pres, .id = "ontogeny") %>%
  filter(pres == 1) %>%
  dplyr::select(tax, obsid, inventory, ontogeny) %>%
  filter(tax %in% species)

selectspecies <- Species %>%
  group_by(tax, inventory, obsid, ontogeny) %>%
  tally() %>% # View()
  filter(n >= obsthreshold) %>% # Filter 2: obsthreshold observations per tax/inventory/ontogeny/obsid
  group_by(tax, inventory, obsid) %>% # This is already the case, but making it excplicit
  filter(any(ontogeny==1) & any(ontogeny == 2)) %$%  # Filter 1: Both ontogenetic stages present within grouping (any obsid)
  tax %>%
  unique()


subsetTax <- function(tx, Data){
  Data %<>% dplyr::filter(tax == tx) %>%
    droplevels() #! Otherwise e.g. all taxa and plots will be carried along.
}

taxtables_small <- lapply(selectspecies, subsetTax, Data = Small_pres)
names(taxtables_small) <- selectspecies

taxtables_big <- lapply(selectspecies, subsetTax, Data = Big_pres)
names(taxtables_big) <- selectspecies

rm(Small_pres, Big_pres)


# Taxonomic treatment -----------------------------------------------------
## For other projects and harmonization across inventories, Pinus uncinata has been lumped into Pinus mugo. But with the current subset only Spain has records of P. mugo, which are all P. uncinata.
## But be aware that tax.id remains the same!
renameSpecies <- function(taxtablelist, from, to){
  # taxon <- str_replace(to, "\\.", " ")
  # tax.genus <- str_extract(to, "^\\w*")
  # tax.spec <- str_extract(to, "\\w*$")

  names(taxtablelist)[names(taxtablelist) == from] <- to
  taxtablelist[[to]]$tax <- to
  return(taxtablelist)
}

taxtables_small %<>% renameSpecies(from = "Pinus.mugo", to = "Pinus.uncinata")
taxtables_big %<>% renameSpecies(from = "Pinus.mugo", to = "Pinus.uncinata")



# Join predictors -----------------------------------------------------------
joinEnv <- function(Data, joinby = "envjoinid"){
  # joinid_Env <- coalesce(Env$plotobsid, Env$plotid, Env$clusterid) ## beautiful but unnecessary
  refsys <- sf::st_crs(Env)
  Data %<>%
    droplevels() %>%
    left_join(dplyr::select(Env, -one_of(setdiff(names(Data), joinby))) %>% mutate_at(joinby, as.character), by = joinby) %>% # This way round because there is NAs in the ids in Env.
    sf::st_as_sf(crs = refsys)
}

taxtables_small %<>% lapply(joinEnv, joinby = "plotid") # !!! joinby = "plotid", not the default, only here, because unique environments per plot are sought (see above).
taxtables_big %<>% lapply(joinEnv, joinby = "plotid") # Warnings about 'Unknown columns' are A-Ok.
## warnings from one_of are cool!

# NA handling -----------------------------------------------------------

sumNAcols <- function(Data){
  colSums(is.na(Data))
}

sumEmptyGeomRows <- function(Data){
  sum(st_is_empty(Data))
}

dropNArows <- function(Data, exceptcol){
  nrow_in <- nrow(Data)
  Wtfcompletecasesisshitty <- Data[, !names(Data) %in% exceptcol]
  if(is(Data, "sf")) Wtfcompletecasesisshitty %<>% st_drop_geometry()
  Data <- Data[complete.cases(Wtfcompletecasesisshitty),]
  cat('Dropping', paste(nrow_in - nrow(Data)), 'of', paste(nrow_in), 'rows.\n')
  return(Data)
}

dropEmptyGeomRows <- function(Data){
  nrow_in <- nrow(Data)
  Data %<>% filter(!st_is_empty(geometry))
  cat('Dropping', paste(nrow_in - nrow(Data)), 'of', paste(nrow_in), 'rows.\n')
  return(Data)
}

## Explore empty Geom rows
# taxtables_small %>% lapply(sumEmptyGeomRows)
# taxtables_big %>% lapply(sumEmptyGeomRows)

## Usually there are just two plotids in DE_BWI which do not correspond to anything, including geometries in Env table.
# taxtables_small %<>% lapply(dropEmptyGeomRows)
# taxtables_big %<>% lapply(dropEmptyGeomRows)

## Explore NA cols
# taxtables_small %>% lapply(sumNAcols)
# taxtables_big %>% lapply(sumNAcols)

# taxtables_big[[1]] %>% dplyr::filter(is.na(wwpi_cop))


## Drop NA rows.
##  Note that they might already be handled above by imputation.
##  Inventory-specific ids are excluded from dropping
exceptcol_small <- c('clusterobsid')
exceptcol_big <- c('clusterobsid', 'plotarea_le40mm', 'plotarea_gr39mmle100m', 'plotarea_gr99mm')

taxtables_small %<>% lapply(dropNArows, exceptcol = exceptcol_small)
taxtables_big %<>% lapply(dropNArows, exceptcol = exceptcol_big)


# Geometry to predictors --------------------------------------------------
## This also drops sf attributes!

convertGeometry <- function(Data){
  LonLat <- st_coordinates(st_geometry(Data$geometry))
  Data %<>%
    mutate(lon = LonLat[,1], lat = LonLat[,2]) %>%
    st_drop_geometry()
}

taxtables_small %<>% lapply(convertGeometry)
taxtables_big %<>% lapply(convertGeometry)


# Subset to obsids where both life stages have been assessed --------------
# This is much faster than dplyr::group_by approaches

intersectID <- function(S, B, idname){
  base::intersect(dplyr::pull(S, idname), dplyr::pull(B, idname))
}

subsetToID <- function(Table, id){
  Table[Table$plotobsid %in% id,]
}

plotobsids_intersect <- mcmapply(intersectID,
                                 taxtables_small,
                                 taxtables_big,
                                 MoreArgs = list(idname = 'plotobsid'),
                                 mc.cores = n_cores)

taxtables_small <- mcmapply(subsetToID, taxtables_small, plotobsids_intersect, SIMPLIFY = F, mc.cores = n_cores)
taxtables_big <- mcmapply(subsetToID, taxtables_big, plotobsids_intersect, SIMPLIFY = F, mc.cores = n_cores)



# Combine ontogenetic stages into one dataset -------------------------------------
## expects two, returns one data.frameoid
combineSmallBig <- function(S, B) {
  S %<>%
    mutate_at(c('methodid'), function(.) paste('small', ., sep = '_'))
  B %<>%
    mutate_at(c('methodid'), function(.) paste('big', ., sep = '_'))
  
  dplyr::bind_rows(small = S, big = B, .id = 'ontogeny')
}

taxtables_comb <- mapply(FUN = combineSmallBig,
                         taxtables_small,
                         taxtables_big,
                         SIMPLIFY = F)

## Shortcut
# saveRDS(taxtables_comb, "Data/taxtables_comb_temp.rds")
# taxtables_comb <- readRDS("Data/taxtables_comb_temp.rds")



# Subset data to reasonable intersections  ----------------------------------
# the data are divided by species anyway!

## Inspect
# ## There are 8 additional small tree observations in DE_BWI and 32411 less small tree observations in some species (De_BWI1987).
# ## Subsetting to intersecting plotobsids is done to keep the data balanced anyway.
# 
# 
# setdiffID <- function(S, B, idname){
#   diff1 <- base::setdiff(dplyr::pull(S, idname), pull(B, idname))
#   diff2 <- base::setdiff(dplyr::pull(B, idname), pull(S, idname))
#   list(length1diff2 = length(diff1), length2diff1 = length(diff2))
# }
# 
# 
# ## Inspect how many diverging plots there are intraspecifically across big and small tables.
# mcmapply(setdiffID,
#          taxtables_small,
#          taxtables_big,
#          MoreArgs = list(idname = 'plotobsid'),
#          mc.cores = n_cores)


filterFirstLastRepeated <- function(D){
  D %<>%
    mutate(obsid_num = as.numeric(as.factor(obsid))) %>% ## Have to be alphabetically ordered as in DE_BWI_2002, DE_BWI_2012
    
    ## hacky way of efficiently getting logicals for last and first obsid
    group_by(inventory) %>%
    mutate(firstobsid = min(obsid_num), lastobsid = max(obsid_num)) %>% 
    mutate(firstobsid = obsid_num == firstobsid, lastobsid = obsid_num == lastobsid) %>%
    ungroup() %>%
    
    dplyr::filter(firstobsid | lastobsid) %>%
    dplyr::select(-obsid_num, -firstobsid, -lastobsid) %>%
    droplevels() %>%
  
    ##
    group_by(inventory, plotid) %>%
    mutate(n_surveys = n_distinct(obsid)) %>%
    ungroup() %>%
    
    filter(n_surveys > 1) %>%
    dplyr::select(-n_surveys) %>%
    droplevels()
  
  return(D)
}


## 1. Data have been subset to all plotids which have been surveyed for both ontogeny levels within plotobsids (regardless of presence/absence)
## !!! On purpose before subsetting to first/last


## 2. Subset to those observations, which were in the first or last species-specific survey per inventory, and which are also repeated within plotid
## !!! On purpose before subsetting to intersection of plotids across obsids. (Subsetting to first last is done for effect comparison).
taxtables_comb_timeintersect <- taxtables_comb %>% mclapply(filterFirstLastRepeated, mc.cores = n_cores)  # for some reason, the plots in two clusters aren't repeated in ES_IFN


# ## base way of doing the repetition intersect, probably faster!
# intersectAcrossObsid <- function(D, idname){
#   ids <- split(dplyr::pull(D, idname), D$obsid)
#   inventoryprefix <- str_extract(names(ids), "^[a-zA-Z]*_[a-zA-Z]*")
#   ids_inv <- split(ids, inventoryprefix) # A list of inv/obsid/c(plotids)
#   unlist(lapply(ids_inv, function(i) Reduce(base::intersect, i)))
# }
# 
# intersectingplotids <- mapply(intersectAcrossObsid,
#                               taxtables_comb,
#                               MoreArgs = list(idname = 'plotid'))
# 
# taxtables_comb_timeintersect <- mapply(subsetToPlotid, taxtables_comb, intersectingplotids, SIMPLIFY = F)


## Check how much of the observations are left
BeforeAfter <- bind_rows(lapply(taxtables_comb, nrow), lapply(taxtables_comb_timeintersect, nrow))


# tabulatePlots <- function(t){
#   t %<>%
#     bind_rows() %>%
#     dplyr::filter(inventory == "ES_IFN")
#   table(t$plotid, t$clusterid)
# }

# N <- tabulatePlots(taxtables_comb) > 0
# all(apply(N, 1, sum) == 1)



# Transform variables -----------------------------------------------------
# transformVariables <- function(Data){
#   dplyr::mutate(Data,
#                 cn_esdacc_sqrt = sqrt(cn_esdacc), ## sqrt()-transform C/N because of inherent skewness.
#                 cn_esdacc_log = log(cn_esdacc),
#                 sand_esdact_log = log(sand_esdact))
# }
# 
# predictor_transformed <- c('cn_esdacc_sqrt', 'cn_esdacc_log', 'sand_esdact_log')
# 
# 
# taxtables_comb %<>% lapply(transformVariables)
# taxtables_comb_timeintersect %<>% lapply(transformVariables)
# Env %<>% transformVariables # Do the same with Env for obtaining center and scale parameters!


# Scaling ------------------------------------------------------------------------
betamin <- .Machine$double.eps
betamax <- 1 - betamin # sprintf("%.16f", betamax)

getScale <- function(Data, colnames){
  Data %<>% mutate_at(colnames, as.numeric)
  Mean <- Data %>% summarize_at(colnames, mean, na.rm = T)
  SD <- Data %>% summarize_at(colnames, sd, na.rm = T)
  Min <- Data %>% summarize_at(colnames, min, na.rm = T)
  Max <- Data %>% summarize_at(colnames, max, na.rm = T)
  bind_rows('mean' = Mean, 'sd' = SD, "min" = Min, "max" = Max, .id = 'scalepar')
}

scaleCols <- function(Data, colnames){
  Data %<>%
    mutate_at(colnames, as.numeric) %>%
    mutate_at(colnames, list("01" = ~ rescale(., to = c(betamin, betamax)))) %>% # Scaling to 0…1 (from range) with 01 suffix
    mutate_at(across(all_of(colnames),
                     .fns = ~ c(scale(.x)))) # Ordinary scaling without suffix
}

scaleColsGrand <- function(Data, colnames, Scaling = Scale_comb){
  # Get "from" parameters for scaling
  min <- unlist(Scaling[3, colnames])
  max <- unlist(Scaling[4, colnames])
  m <-  unlist(Scaling[1, colnames])
  s <- unlist(Scaling[2, colnames])
  
  Data %<>%
    mutate_at(colnames, as.numeric) %>%
    mutate(across(.cols = all_of(colnames),
                  .fns = ~ rescale(.x, to = c(betamin, betamax), from = c(min[cur_column()], max[cur_column()])),
                  .names = "{col}_01")) %>% # Scaling to 0…1 with 01 suffix
    mutate(across(all_of(colnames),
                  .fns = ~ c(scale(.x, center = m[cur_column()], scale = s[cur_column()])))) # Ordinary scaling without suffix
  return(Data)
}


# ### Special scaling to fixed ranges for time.
# getTimeExtreme <- function(tables, fun = "min", selectinventory = NULL){
#   if (!is.null(selectinventory)) tables <- lapply(tables, function(t)
#     dplyr::filter(t, inventory == selectinventory))
#   get(fun)(sapply(tables, function(X)  get(fun)(X$time, na.rm = T)), na.rm = T)
# }


getTimeRanges_ES <- function(tables){
  tables %>%
    bind_rows() %>%
    dplyr::filter(inventory == "ES_IFN") %>%
    group_by(clusterid) %>%
    summarize(min = min(time, na.rm = T), max = max(time, na.rm = T), surveys = n_distinct(obsid, na.rm = T), n = n()) %>%
    mutate(period = as.numeric(max) - as.numeric(min))  %>%
    ungroup()
}


scaleTime <- function(tables, refperiodyears = 25){
  
  refperioddays <- refperiodyears * 365.25
  refperiodrange <- round(c(-refperioddays*0.5, refperioddays*0.5))
  
  tables %<>%
    bind_rows(.id = "splitid") %>%
    mutate(time_unscaled = time, time = as.numeric(time)) %>%
    
    group_by(inventory, tax) %>%
    mutate(invstart = min(time), invend = max(time), invtimerange = invend - invstart, invtimecenter = (invstart + invend)*0.5) %>%
    ungroup() %>%
    
    group_by(clusterid, tax) %>%
    mutate(clusterstart = min(time), clusterend = max(time), clustertimerange = clusterend - clusterstart, clustertimecenter = (clusterstart + clusterend)*0.5) %>%
    ungroup() %>%
    
    ## 1. center time per the level, that the time slopes are fitted inventory/cluster(for ES)/taxon (we are only interested in period, not actual time).
    mutate(time = if_else(inventory == "ES_IFN", time - clustertimecenter, time - invtimecenter)) %>%
    
    ## 2. scale time per ref period
    ## time scaled  on -0.5…0.5 from reference time range ...
    ## range is chosen, so that the model effect size (at time increase of 1) corresponds to the actual reference time range.
    mutate(time = rescale(time, to = c(-0.5, 0.5), from = refperiodrange)) %>%
    split(f = .$splitid)
}


#### TIME SCALING

## Periods in Spain are very unbalanced across geographical clusters.
# TimeRange_ES <- getTimeRange_ES(taxtables_comb_timeintersect)
# maxperiod_ES <- TimeRange_ES[match(max(TimeRange_ES$period), TimeRange_ES$period),]

## This is why the time is centered per plot and scaled by maximum plot time range:
taxtables_comb %<>% scaleTime()
taxtables_comb_timeintersect %<>% scaleTime()


## The maximal time extent per plot is in Spain!
# SE_NFI: 3653 / 365 = 10yrs
# DE_BWI: 9132 / 365 = 25yrs
# ES_IFN: 12418 / 365 = 34yrs 


## in overall periods
## ES IFN: 1980-2017
## DE_BWI: 1987-2012 (25); all equidistant
## SE_NFI: 2003-2017 (14)

## old method:
# refcountry <- "ES_IFN" ## ... for one reference time period.
# timemin_ref <- getTimeExtreme(taxtables_comb, "min", refcountry)
# timemax_ref <- getTimeExtreme(taxtables_comb, "max", refcountry)



### Save data with only time scaled.
saveRDS(taxtables_comb, 'Data/taxtables_comb_unscaled.rds')
saveRDS(taxtables_comb_timeintersect, 'Data/taxtables_comb_timeintersect_unscaled.rds')
# saveRDS(taxtables_comb, 'Data/taxtables_comb_unscaled_v2.rds', version = 2)



scalecols <- c(predictor_select) # , predictor_transformed) "lon", "lat", only in taxtables, not in env

### Scale per table:
## Save scaling parameters for per table scaling!
# scales_comb <- lapply(taxtables_comb, getScale, colnames = scalecols)
# saveRDS(scales_comb, 'Data/Scalings_species_comb.rds')
## Scale
# taxtables_comb %<>% lapply(scaleCols, colnames = scalecols)

### Scale per grand mean/sd:
## in order to have all trees (all species, small, big; big/small combined scaling would for the most part be the case anyway)
## and all different sampling methods/protocols/scope on the same scale.

## Env is already unique per plotid!
# Env_unique <- Env[any(duplicated(Env$plotid)),]


Scale_comb <- getScale(as.data.frame(Env), colnames = scalecols)

saveRDS(Scale_comb, 'Data/Scaling_grand_comb.rds')


## Scale scalecols to grand parameters.
scalecols_tax <- scalecols
taxtables_comb %<>% lapply(scaleColsGrand, colnames = scalecols_tax)
taxtables_comb_timeintersect %<>% lapply(scaleColsGrand, colnames = scalecols_tax)



## Split data
taxtables_split <- lapply(taxtables_comb, function(taxtable) split(taxtable, taxtable$ontogeny))
taxtables_small <- lapply(taxtables_split, function(t) t$small)
taxtables_big <- lapply(taxtables_split, function(t) t$big)
rm(taxtables_split)

taxtables_split_timeintersect <- lapply(taxtables_comb_timeintersect, function(taxtable) split(taxtable, taxtable$ontogeny))
taxtables_small_timeintersect <- lapply(taxtables_split_timeintersect, function(t) t$small)
taxtables_big_timeintersect <- lapply(taxtables_split_timeintersect, function(t) t$big)
rm(taxtables_split_timeintersect)

# Save taxtables ---------------------------------------------------------------
saveRDS(taxtables_comb, 'Data/taxtables_comb.rds')
saveRDS(taxtables_small, 'Data/taxtables_small.rds')
saveRDS(taxtables_big, 'Data/taxtables_big.rds')

saveRDS(taxtables_comb_timeintersect, 'Data/taxtables_comb_timeintersect.rds')
saveRDS(taxtables_small_timeintersect, 'Data/taxtables_small_timeintersect.rds')
saveRDS(taxtables_big_timeintersect, 'Data/taxtables_big_timeintersect.rds')



# Generate taxtables with environment at presences only -------------------

obsthreshold <- 20 # minimal n_obs within species/inventory/obsid/ontogeny

## Shortcut:
# taxtables_comb_timeintersect  <- readRDS('Data/taxtables_comb_timeintersect.rds')

subsetToPresence <- function(Data){
  Data %>% 
    dplyr::filter(as.logical(pres)) %>%
    droplevels()
}

taxtables_pres <- taxtables_comb_timeintersect %>% lapply(subsetToPresence)

## Filter species again, after all the subsetting.
purgeInventories <- function(Data){
  
  ## Character of inventories, where the threshold is not met in one obsid.
  inventory_subthreshold <- group_by(Data, inventory, obsid, ontogeny) %>%
    tally %>%
    dplyr::filter(n < obsthreshold) %$% inventory %>%
    unique()
  ## Character giving the inventories, where there is only one obsid.
  inventory_single <- dplyr::group_by(Data, inventory, obsid) %>%
    tally() %>% tally(wt = n()) %>% filter(n < 2) %$% inventory
  
  filter(Data, !(inventory %in% c(inventory_subthreshold, inventory_single)))
}

taxtables_pres %<>% lapply(purgeInventories) %>% .[sapply(., function (x) nrow(x) > 0)]

#### Save
saveRDS(taxtables_pres, "Data/taxtables_pres_thresholdsubset.rds")

#### Anonymize data taxtables_pres
## 1. Get rid of coordinates
taxtables_pres %<>% lapply(function(Tab) dplyr::select(Tab, -lon, -lat))
## 2. Synonymize plotids, clusterids
synonymizeTab <- function(Tab) {
  Tab %<>%
    select(-envjoinid, -plotobsid, -clusterobsid) %>%
    group_by(clusterid) %>%
    mutate(clusterid = if_else(inventory != "ES_IFN", paste(inventory, cur_group_id(), sep = "_"), clusterid)) %>%
    group_by(plotid) %>%
    mutate(plotid = paste(inventory, clusterid, cur_group_id(), sep = "_")) %>%
    ungroup()
  return(Tab)
}

taxtables_pres %<>% lapply(synonymizeTab)
  

saveRDS(taxtables_pres, "Data/taxtables_pres_thresholdsubset_anon.rds")

## Check again.
# taxtables_pres %>% lapply(function(D) table(D$obsid, D$ontogeny))


