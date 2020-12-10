# The Library ---------------------------------------------------------
## Wrangling.
library(tidyverse) # ggplot2, purrr, tibble, dplyr, tidyr, stringr, readr, forcats
library(magrittr)
library(glue)
library(stringi)

## Fitting.
library(glmmTMB) # devtools::install_github("glmmTMB/glmmTMB/glmmTMB")
library(parallel)

## Check and predict.
library(ggeffects)
library(effects)
library(DHARMa)
library(qgam)
library(psych)

# Tabulating.
library(sjPlot)
library(emmeans)
library(lsmeans)
library(multcomp)
library(extrafont)
library(colorspace)
library(lubridate)

## Plotting.
library(ggplot2)
library(cowplot)
library(gridGraphics)

runid <- glue("{stri_rand_strings(1, 5)}") # runid <- "DEF"
runid
set.seed(1)
n_cores <- 3
writeLines(capture.output(sessionInfo()), glue("Session/{runid}_sessionInfo__{format(Sys.time(), '%Y-%m-%d')}.txt"))
saveRDS(sessionInfo(), glue("Session/{runid}_sessionInfo__{format(Sys.time(), '%Y-%m-%d')}.rds"))


# Set working directory ---------------------------------------------------
setwd(rprojroot::find_root(rprojroot::criteria$is_git_root))

# Source ------------------------------------------------------------------
## Presumed to have run:
# source('Prepare data.R')

## The objects in this script are presumed to be loaded
source('Theme.R') # Load ggplot themes, colors.

# Load data ---------------------------------------------------------------

taxtables_pres <- readRDS('Data/taxtables_pres_thresholdsubset_anon.rds')


# Requirements for plotting -------------------------------------------------
objectify <- function(plotcall){
  par(xpd = NA, # switch off clipping, necessary to always see axis labels
      bg = "transparent") # switch off background to avoid obscuring adjacent plot
  plotcall
  return(recordPlot())
}


# Never stop wrangling ----------------------------------------------------

## Split by inventory for inventory level models.
splitByInventory <- function(Data){
  Data %>%
    split(.$inventory)
}

invtables <- taxtables_pres %>% lapply(splitByInventory) %>% unlist(recursive = F)



## Select final niche dimensions.
nichedim <- c('gdd0_aclim_01',
              'cwbYear_aclim_01')


#############################################################################
# ANALYSIS ------------------------------------------------------------------
#############################################################################

# Fit the model ----------------------------------------------------------------
modelname <- "glmmTMB_beta_ontogeny_time"
modelfixef <- "time*ontogeny"
modeldispef <- "time*ontogeny"

fitL <- function(taxname, envname, invname, glmmTMBfamily = "beta_family",
                 fixef = modelfixef, extraeffect, dispef = modeldispef) {
  glmmTMBfamily <- glmmTMBfamily
  rhs <- glue("{fixef} + {extraeffect}")
  tablename <-  glue("{taxname}.{invname}")
  Data <- invtables[[tablename]]
  disppercluster <- F # invname == "ES_IFN"

  formula <- glue("{envname} ~ {rhs}") %>% as.formula()
  dispformula <- if(disppercluster) as.formula(glue("~ {dispef}*clusterid")) else as.formula(glue("~ {dispef}"))

  fit <- tryCatch(
    {
      f <- glmmTMB::glmmTMB(formula,
                      data = Data,
                      family = glmmTMBfamily,
                      dispformula = dispformula,
                      REML = F, # !
                      control = glmmTMBControl(optCtrl=list(iter.max=1e5, eval.max=1e5)))
    },
    # warning = function(w) message(w), ## don't handle warnings
    error = function(e) return(e),
    finally = {}
    )
  
  ## just add attributes to the fit
  fit$inventory <- invname
  fit$env <- envname
  fit$tax <- taxname
  
  return(fit)
}


C <- Comb <- expand.grid(datacomb = names(invtables), env = nichedim, stringsAsFactors = F) %>%
  mutate(tax = str_extract(datacomb, "[^.]*\\.[^.]*"), inv = str_extract(datacomb, "([^\\.]+$)")) %>%
  mutate(modelextraeffect = ifelse(inv == "ES_IFN", "clusterid", "1")) %>% # "(1 | clusterid)"
  dplyr::arrange(inv, env, tax)

saveRDS(Comb, glue("Data/{runid}_Combinations_{modelname}.rds"))

fits <- mcmapply(fitL, Comb[["tax"]], Comb[["env"]], Comb[["inv"]], extraeffect = Comb[["modelextraeffect"]],
                 SIMPLIFY = F, mc.cores = n_cores) %>% 
  setNames(paste(Comb[["tax"]], Comb[["env"]],  Comb[["inv"]], sep = "__"))

# fits %>% lapply(summary)



# Save and load fits ------------------------------------------------------
saveRDS(fits, glue("Fits/{runid}_Fit_{modelname}_{format(Sys.time(), '%Y-%m-%d')}.rds"))



# runid <- "DEF"
loadrunid <- runid

loadrundate <- format(Sys.time(), '%Y-%m-%d')
# loadrundate <- "2020-11-02" # format(Sys.time(), '%Y-%m-%d')

modelname <- "glmmTMB_beta_ontogeny_time"
modelfixef <- "time*ontogeny"
 
# fits <- readRDS(glue("Fits/{loadrunid}_Fit_{modelname}_{loadrundate}.rds"))
# Comb <- readRDS(glue("Data/Combinations_{modelname}.rds"))


# Subset to converged models ---------------------------------
isglmmTMB <- fits %>%  sapply(is, class2 = "glmmTMB") %>% as.logical() %>% replace_na(F)
isconverged <- fits %>% sapply(function(m) m$sdr$pdHess == T) %>% as.logical() %>% replace_na(F) # The equation hack always returns something.
all(isconverged)

selectmodel <- isglmmTMB & isconverged
Drop <- Comb[!selectmodel,]
Drop
C <- Comb[selectmodel,]
table(selectmodel)
fits <- fits[selectmodel]

# In DEF none dropped!

#############################################################################
# Posterior checks ---------------------------------------------------------
#############################################################################

#### Solved problems with normal regression residuals:
## Former problem with a ontogeneny*time*inventory model (1 | inventory/clusterid): blocked residuals per inventory, probably due to regularization in random effects. 
## Splitting data by inventory and fitting ontogeny*time+(1|clusterid) solves this for DE and SE.
## But particularly in ES, there still seem to be two groups with different dispersal (blocks).
## Groups in residuals are basically interaction(ontogeny, time), this is why 4 in DE_BWI.

#### ES problem:
## There was some bias of residuals/response against time because in Spain different geographical clusters have been surveyed across differen periods, e.g. 1. mountaineuos region 1980-1995, 2. coastal region 1990-2015.
## This has been dealt with centering clusters within inventory.



plotCluster <- function(d) {
  D <- d$fittedModel$frame
  boxplot(d$scaledResiduals ~ D$clusterid)
}

# Residual checks --------------------------------------------------------
dharmas <- fits %>% mclapply(simulateResiduals, mc.cores = n_cores)
saveRDS(dharmas, glue("Residuals/{runid}_Residuals_{modelname}_{format(Sys.time(), '%Y-%m-%d')}.rds"))
# dharmas <- readRDS("Fits/{loadrunid}_Residuals_{modelname}_{loadrundate}.rds"))

dharmaplots <- dharmas %>%
  lapply(function(d) objectify(plot(d, quantreg = F)))

clusterdharmaplots <- lapply(1:length(dharmas),
                             function(i) recalculateResiduals(dharmas[[i]], group = invtables[[C[i,"tax"]]]$clusterid)) %>%
  lapply(function(d) objectify(plot(d, quantreg = F)))

# dharmaplots_clusterbox <- dharmas %>%
#   lapply(function(d) objectify(plotCluster(d))) # somehow plotResiduals(d, form = ~ clusterid) complains.

n_nichedim <- length(unique(C[["env"]]))
n_species <- length(unique(C[["tax"]]))
dharmalabels <- sapply(dharmas, function(d) paste(d$fittedModel$tax, d$fittedModel$env, d$fittedModel$inv))

dharmaplotgrid <- cowplot::plot_grid(plotlist = dharmaplots,
                            # align = "hv",
                            labels = dharmalabels,
                            axis = "bl",
                            ncol = n_nichedim, # No. of selected nichedims.
                            greedy = F)

# clusterdharmaplotgrid <- cowplot::plot_grid(plotlist = clusterdharmaplots,
#                                      # align = "hv",
#                                      labels = dharmalabels,
#                                      axis = "bl",
#                                      ncol = n_nichedim, # No. of selected nichedims.
#                                      greedy = F)

# clusterdharmaboxgrid <- cowplot::plot_grid(plotlist = dharmaplots_clusterbox,
#                                             # align = "hv",
#                                             labels = dharmalabels,
#                                             axis = "bl",
#                                             ncol = n_nichedim, # No. of selected nichedims.
#                                             greedy = F)

cowplot::save_plot(glue("Residuals/{runid}_Residuals_{modelname}_{format(Sys.time(), '%Y-%m-%d')}.pdf"),
          dharmaplotgrid,
          device = "pdf",
          base_width = 160 * n_nichedim,
          base_height = 120 * n_species,
          units = "mm", limitsize = F)

# cowplot::save_plot(glue("Residuals/{runid}_Residuals_clusterid_{modelname}_{format(Sys.time(), '%Y-%m-%d')}.pdf"),
#                    clusterdharmaplotgrid,
#                    device = "pdf",
#                    base_width = 160 * n_nichedim,
#                    base_height = 120 * n_species,
#                    units = "mm", limitsize = F)

# cowplot::save_plot(glue("Residuals/{runid}_Residuals_clusterboxes_{modelname}_{format(Sys.time(), '%Y-%m-%d')}.pdf"),
#                    clusterdharmaboxgrid,
#                    device = "pdf",
#                    base_width = 160 * n_nichedim,
#                    base_height = 120 * n_species,
#                   units = "mm", limitsize = F)



# Drop models based on residuals ------------------------------------------

droppedmodel <- c("Alnus.incana__gdd0_aclim_01__DE_BWI", "Populus.tremula__gdd0_aclim_01__SE_NFI", "Salix.caprea__gdd0_aclim_01__SE_NFI")

dropmodelfromfits <- names(fits) %in% droppedmodel
fits <- fits[!dropmodelfromfits]



#############################################################################
# ESTIMATES ----------------------------------------------------------------
#############################################################################


# Effects and inference -----------------------------------------------------------

## Notes on glmmTMB and emmeans:
## 1. The dispersion model has a log link, effects are on a log scale which is fine!
##    Nevertheless the same tests are applied. This way contrasts can be tested against 0, and diminishing variance is indicated by beta < 0. Difference contrasts will turn into ratio contrasts.  
##    Alternative tried before was:
##     — type = "response", but only AFTER contrasts from emmeans.
##     - This does not seem to work with predictions in emmtrends. But a link function can be created to do it manually.
## 2. The dispersion parameters correspond to variance = sd^2. Consider sqrt().
## 3. Dispersion (here variance) model can be extracted from fits with emmeans(fit, component = "disp") as in https://github.com/glmmTMB/glmmTMB/blob/78856c092e286f4d449af707bd7fdfaf638957c3/glmmTMB/R/emmeans.R
## 4. As is, if average over inventories are desired (averageoverinventories == TRUE), they are weighted proportionally to n observations (emmeans(…, weights ="proportional"))


## Create link functions for transforming emmGrids (including effects, CI).
## Absolute for effect magnitudes.
abslink <- make.link("identity")
abslink$linkfun <- function(mu) abs(mu)
abslink$name <- "abs"


## Returns a data.frame of ...
### 1. the ontogeny effect:
## (small intercept - big intercept) [hence: revpairwise!] by inventory,
## The null whether == 0 is tested (disp model on log scale, which means on response scale == 1).
### 2. time slopes (time effects) by ontogeny and inventory,
## The null whether == 0 is tested (disp model on log scale, which means on target scale == 1) .

getEffects <- function(fit,
                       glmmTMBcomponent = "cond", # Which model component to use, dispersion or "conditional" (µ)? Passed on to https://github.com/glmmTMB/glmmTMB/blob/78856c092e286f4d449af707bd7fdfaf638957c3/glmmTMB/R/emmeans.R Alternatives: "cond", "disp".
                       attime = NULL, # c(-1, 0, 1), # for the ontogeny effect
                       ontogenylevel = c("small", "big"), # for the time effects
                       abs = FALSE,
                       transform = "response", # alternative: "none"
                       averageoverinventories = FALSE) {
  
  glmmTMBcomponent <- glmmTMBcomponent ## https://stackoverflow.com/questions/18136720/why-missing-and-default-arguments-are-not-working-in-functions-called-by-lapp
  
  
  oneinventory <- is.null(fit$frame$inventory)
  if (oneinventory & is.null(attime)) {

    ont <- emmeans::emmeans(fit,
                            specs = revpairwise ~ ontogeny,
                            transform = transform,
                            component = glmmTMBcomponent)$contrasts
    timetrend <- emmeans::emtrends(fit,
                                   specs = ~ ontogeny,
                                   var = "time",
                                   transform = transform,
                                   component = glmmTMBcomponent)
    
    timetest <- emmeans::test(timetrend,
                              null = 0,
                              joint = F,
                              by = c("ontogeny")) # "by" ensures independent tests!
    
  } else if(oneinventory) {
    ont <- emmeans::emmeans(fit,
                            specs = revpairwise ~ ontogeny | time,
                            at = list(time = attime),
                            transform = transform, 
                            component = glmmTMBcomponent)$contrasts
    timetrend <- emmeans::emtrends(fit,
                                   specs = ~ ontogeny,
                                   var = "time",
                                   transform = transform,
                                   component = glmmTMBcomponent)
    
    timetest <- emmeans::test(timetrend,
                              null = 0,
                              joint = F,
                              by = c("ontogeny")) # "by" ensures independent tests!
  } else if (averageoverinventories){
    ont <- emmeans::emmeans(fit,
                            specs = revpairwise ~ ontogeny | time,
                            at = list(time = attime),
                            transform = transform,
                            weights = "proportional", # For average predictions, here over inventories. Alternatives: here "outer" == "proportional", "equal"
                            component = glmmTMBcomponent
                            )$contrasts
    timetrend <- emmeans::emtrends(fit,
                                   specs = ~ ontogeny,
                                   var = "time",
                                   weights = "proportional", # For average predictions, here over inventories. Alternatives: here "outer" == "proportional", "equal"
                                   component = glmmTMBcomponent)
    
    timetest <- emmeans::test(timetrend,
                              null = 0,
                              joint = F,
                              by = c("ontogeny")) # "by" ensures independent tests!
  } else {
    ont <- emmeans::emmeans(fit,
                            specs = revpairwise ~ ontogeny | inventory:time,
                            at = list(time = attime),
                            transform = transform,
                            component = glmmTMBcomponent
                            )$contrasts
    timetrend <- emmeans::emtrends(fit,
                                   specs = ~ ontogeny | inventory,
                                   var = "time",
                                   transform = transform,
                                   component = glmmTMBcomponent)
    
    timetest <- emmeans::test(timetrend,
                              null = 0,
                              joint = F,
                              by = c("inventory", "ontogeny")) # "by" ensures independent tests!
  }
  
  #### Ontogeny effect quantities
  CL_ont <- as.data.frame(confint(ont))[c("lower.CL", "upper.CL")] # level 0.95;
  p_rejneg_ont <- test(ont, side = ">")$p.value
  p_rejpos_ont <- test(ont, side = "<")$p.value
  
  Ont <- cbind(as.data.frame(ont), CL_ont,  p_rejneg = p_rejneg_ont,  p_rejpos = p_rejpos_ont, effect = "ontogeny")
  
  
  #### Ontogeny effect quantities
  CL_time <- as.data.frame(timetrend)[c("lower.CL", "upper.CL")] # level 0.95
  p_rejneg_time <- test(timetrend, side = ">")$p.value
  p_rejpos_time <- test(timetrend, side = "<")$p.value
  
  Time <- cbind(as.data.frame(timetest), CL_time, p_rejneg = p_rejneg_time, p_rejpos = p_rejpos_time, effect = "time") %>%
    dplyr::filter(ontogeny %in% ontogenylevel) %>%
    dplyr::rename(estimate = "time.trend")
  
  
  #### Bind and process both
  Out <- bind_rows(Ont, Time) %>% 
    mutate(signif = p.value < 0.05, signif_rejneg = p_rejneg < 0.05, signif_rejpos = p_rejpos < 0.05) %>%
    mutate(star = symnum(p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")))

  if (abs & !glmmTMBcomponent == "disp") {
    inverted <- Out$estimate < 0
    lower <- Out$lower.CL
    Out$lower.CL[inverted] <- Out$upper.CL[inverted]
    Out$upper.CL[inverted] <- lower[inverted]
    Out$estimate <- abs(Out$estimate)
  }
  
  if(oneinventory) Out$inventory <- fit$inventory
  Out$env <- fit$env
  Out$tax <- fit$tax
  
  return(Out)
}

## Wrappers for legacy reasons:
getOntogenyEffect <- function(fit, ...) dplyr::filter(getEffects(fit, ...), effect == "ontogeny")
getTimeEffect <- function(fit, ...) dplyr::filter(getEffects(fit, ...), effect == "time")


## Tests for difference between effects with emmeans.
## Returns: data.frame
## Null: Equal means of ontogeny (i.e. big-small) and time effect, i.e. time - ontogeny == 0.
## THIS TEST IS DONE ON THE MODEL SCALE, AS IS THE RESULTING DIFF!

testEffectDiff <- function(fit,
                           glmmTMBcomponent = "cond", # Which model component to use, dispersion or "conditional" (µ)? Passed on to https://github.com/glmmTMB/glmmTMB/blob/78856c092e286f4d449af707bd7fdfaf638957c3/glmmTMB/R/emmeans.R Alternatives: "cond", "disp".
                           attime = NULL, # 0,
                           magnitude = FALSE, # If magnitude == TRUE, will test for difference of absolute effects.
                           averageoverinventories = FALSE) {
  
  glmmTMBcomponent <- glmmTMBcomponent ## https://stackoverflow.com/questions/18136720/why-missing-and-default-arguments-are-not-working-in-functions-called-by-lapp
  
  
  oneinventory <- is.null(fit$frame$inventory)
  n_inventories <- if (oneinventory) 1 else length(unique(fit$frame$inventory))
  

  ## Extract timeslope for "small" and 'ontogenyslope', depending on whether to average over inventories.
  if (oneinventory & is.null(attime)) {
    
    timeslope <- emmeans::emtrends(fit,
                                   specs = ~ ontogeny,
                                   var = c("time"),
                                   component = glmmTMBcomponent
    )
    smalltimeslope <- timeslope[2] #! small!
    
    
    ontogenyslope <- emmeans::emmeans(fit,
                            specs = revpairwise ~ ontogeny,
                            component = glmmTMBcomponent)$contrasts
    
  } else if(oneinventory) {
  
    timeslope <- emmeans::emtrends(fit,
                                   specs = ~ ontogeny,
                                   var = c("time"),
                                   component = glmmTMBcomponent
    )
    smalltimeslope <- timeslope[2] #! small!
    
    ontogenyslope <- emmeans::emmeans(fit,
                                      specs = revpairwise ~ ontogeny | time,
                                      at = list(time = attime),
                                      component = glmmTMBcomponent
                                      )$contrast
    
  } else if (averageoverinventories) {
    timeslope <- emmeans::emtrends(fit,
                                   specs = ~ ontogeny,
                                   var = c("time"),
                                   weights = "proportional", # For average predictions, here over inventories. Alternatives: here "outer" == "proportional", "equal"
                                   component = glmmTMBcomponent
                                   )
    
    smalltimeslope <- timeslope[2] #! small!
    
    ontogenyslope <- emmeans::emmeans(fit,
                                      specs = revpairwise ~ ontogeny | time,
                                      at = list(time = attime),
                                      weights = "proportional", # For average predictions, here over inventories. Alternatives: here "outer" == "proportional", "equal"
                                      component = glmmTMBcomponent #, # transform = "response"
                                      )$contrast

  } else {
    
    timeslope <- emmeans::emtrends(fit,
                                   specs = ~ ontogeny | inventory,
                                   var = c("time"),
                                   component = glmmTMBcomponent #, transform = "response"
                                   )
    smalltimeslope <- timeslope[(1:n_inventories) * 2] #! small!
    
    ontogenyslope <- emmeans::emmeans(fit,
                                      specs = revpairwise ~ ontogeny | inventory:time,
                                      at = list(time = attime),
                                      component = glmmTMBcomponent #, transform = "response"
                                      )$contrast
  }
  
  ## Bind the effect emgrids. This is still on the log scale, even though emmgrid, does not seem to know about it.
  effects <- rbind(smalltimeslope, ontogenyslope, adjust = "none") # Only the difference is tested later, thus no correction necessary here.
  
  ## NOTE: Joint testing would work, but test ignores the side = "<" option :(
  # emmeans::test(effects, null = 0, joint = T)
  
  
  if (magnitude) effects <- regrid(effects, transform = abslink) # This transforms everything adequately through an abs link.
  
  ## Pairwise difference contrasts on the given scale.
  ## Although emmeans complains with a note for the disp model, this is the log-scale!
  diffby <- if (oneinventory) NULL else "inventory"
  diff <- emmeans::contrast(effects,
                            method = "pairwise",
                            by = diffby,
                            adjust = "none")
  
  
  
  
  Diff <- diff %>%
    as.data.frame() %>%
    dplyr::select(-contrast) %>%
    mutate(signif = p.value < 0.05) %>%
    mutate(star = symnum(p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " "))) %>%
    cbind(effect = "diff_time.ontogeny")
  
  if(oneinventory) Diff$inventory <- fit$inventory
  Diff$env <- fit$env
  Diff$tax <- fit$tax
  
  return(Diff)
}


#########################################
#### (I) Get effects on µ.
#########################################

E <- fits %>%
  lapply(getEffects, abs = F,  attime = NULL, ontogenylevel = c("small"), averageoverinventories = FALSE) %>%
  bind_rows(.id = "id") %>%
  dplyr::mutate(env = fct_relevel(env, nichedim)) # For consistent order.

Diff <- fits %>%
  lapply(testEffectDiff) %>% # attime = NULL #!!!
  bind_rows(.id = "id")

Effects <- left_join(E, dplyr::select(Diff, p.value_diff = p.value, signif_diff = signif, id), by = "id")

saveRDS(Effects, glue("Fits/{loadrunid}_Effects_{modelname}_{loadrundate}.rds"))


####################################################
#### (II) Get effects on dispersion/scale
####################################################

# E_disp <- fits %>%
#   lapply(getEffects, glmmTMBcomponent = "disp", abs = F,  attime = NULL, ontogenylevel = c("small"), averageoverinventories = FALSE) %>%
#   bind_rows(.id = "id") %>%
#   dplyr::mutate(env = fct_relevel(env, nichedim)) # For consistent order.
# 
# Diff_disp <- fits %>%
#   lapply(testEffectDiff, glmmTMBcomponent = "disp") %>% # # attime = NULL #!!!
#   bind_rows(.id = "id")
# 
# Diff_magnitude_disp <- fits %>% ## regrid Notes are ok.
#   lapply(testEffectDiff, magnitude = T, glmmTMBcomponent = "disp") %>% # attime = NULL #!!!
#   bind_rows(.id = "id")
# 
# saveRDS(E_disp, glue("Fits/{loadrunid}_Effects_disp_{modelname}_{loadrundate}.rds"))

######################################################################################
# Regression table -----------------------------------------------------------------##
######################################################################################

getData <- function(fit){
  Tax <- dplyr::filter(taxtables_pres[[fit$tax]], inventory == fit$inventory) # !!! use the not yet presence-subsetted plots 
  n_plots <- length(unique(Tax$plotid))
  P <- Tax %>%
    mutate(time_num = as.numeric(time_unscaled)) %>%
    group_by(obsid) %>%
    dplyr::summarize(time_avg = mean(time_num),
                     .groups = 'drop')
  
  periodbetweenaverages <- round((P$time_avg[as.integer(as.factor(P$obsid)) == 2] - P$time_avg[as.integer(as.factor(P$obsid)) == 1])
                                 /365.2425)
  
  
  maxperiod <- round((max(as.numeric(Tax$time_unscaled)) - min(as.numeric(Tax$time_unscaled)))/365.2425)
  year1 <- year(min(Tax$time_unscaled))
  year2 <- year(max(Tax$time_unscaled))
  years <- glue("{year1}–{year2}")
  
  return(data.frame(yr_betweenAvg = periodbetweenaverages, yr_max = maxperiod, n_plots = n_plots, year1 = year1, year2 = year2, years = years))
}

getTitleColumn <- function(fit){
  envname <- recode(fit$env, "gdd0_aclim_01" = "heat sum", "cwbYear_aclim_01" = "water availability")
  invname <- recode(fit$inventory, "SE_NFI" = "Sweden", "DE_BWI" = "Germany", "ES_IFN" = "Spain")
  taxname <- str_replace(fit$tax, fixed("."), " ")
  id <- glue("{fit$tax}__{fit$env}__{fit$inventory}")
  Titles <- data.frame(id = id, env = envname, inventory = invname, tax = taxname)
  return(Titles)
}


printModelRows <- function(fuu) {
  s <- summary(fuu)
  
  Titles <- getTitleColumn(fuu) # One column, will be recycled
  
  if (fuu$inventory != "ES_IFN")  { # 
    C_cond <- t(s$coefficients$cond[c(1:3, nrow(s$coefficients$cond)),c(1,2,4)]) %>%
      as.data.frame() %>%
      dplyr::rename_with(~ paste0(.x, "_cond"))
    
    C_disp <- t(s$coefficients$disp[,c(1,2,4)]) %>%
      as.data.frame() %>%
      dplyr::rename_with(~ paste0(.x, "_disp"))
    
  } else {
    # ontogeny estimate adult
    o_cond_adult <- emmeans::emmeans(fuu, specs = identity ~ ontogeny,
                                     transform = "none", #
                                     component = "cond")$contrasts[1,] %>%
      as.data.frame()
    
    # time estimate adult
    t_cond_adult <- emmeans::emtrends(fuu,
                                      specs = identity ~ ontogeny,
                                      var = "time",
                                      transform = "none",
                                      component = "cond")$contrasts[1,] %>%
      as.data.frame()
    
    
    # ontogeny contrast adult - juvenile, k is the last level
    o_cond_juv <- emmeans::emmeans(fuu, specs = trt.vs.ctrl1 ~ ontogeny,
                                   
                                   transform = "none", #
                                   component = "cond")$contrasts[1,] %>%
      as.data.frame()
    
    # time contrast adult - juvenile
    t_cond_juv <- emmeans::emtrends(fuu,
                                    specs = trt.vs.ctrl1 ~ ontogeny,
                                    
                                    var = "time",
                                    transform = "none",
                                    component = "cond")$contrasts[1,] %>%
      as.data.frame()
    
    C_cond <- data.frame(t(o_cond_adult[,c("estimate", "SE", "p.value")]),
                         t(t_cond_adult[,c("estimate", "SE", "p.value")]),
                         t(o_cond_juv[,c("estimate", "SE", "p.value")]),
                         t(t_cond_juv[,c("estimate", "SE", "p.value")])) %>%
      set_rownames(c("estimate", "se", "p")) %>%
      setNames(c("(Intercept)_cond", "time_cond", "ontogenysmall_cond", "time:ontogenysmall_cond"))
    
    
    
    ##### disp
    # ontogeny estimate adult
    o_disp_adult <- emmeans::emmeans(fuu, specs = identity ~ ontogeny,
                                     transform = "none", #
                                     component = "disp")$contrasts[1,] %>%
      as.data.frame()
    
    # time estimate adult
    t_disp_adult <- emmeans::emtrends(fuu,
                                      specs = identity ~ ontogeny,
                                      var = "time",
                                      transform = "none",
                                      component = "disp")$contrasts[1,] %>%
      as.data.frame()
    
    
    # ontogeny contrast adult - juvenile
    o_disp_juv <- emmeans::emmeans(fuu, specs = trt.vs.ctrl1 ~ ontogeny,
                                   transform = "none", #
                                   component = "disp")$contrasts[1,] %>%
      as.data.frame()
    
    # time contrast adult - juvenile
    t_disp_juv <- emmeans::emtrends(fuu,
                                    specs = trt.vs.ctrl1 ~ ontogeny,
                                    var = "time",
                                    transform = "none",
                                    component = "disp")$contrasts[1,] %>%
      as.data.frame()
    
    C_disp <- data.frame(t(o_disp_adult[,c("estimate", "SE", "p.value")]),
                         t(t_disp_adult[,c("estimate", "SE", "p.value")]),
                         t(o_disp_juv[,c("estimate", "SE", "p.value")]),
                         t(t_disp_juv[,c("estimate", "SE", "p.value")])) %>%
      set_rownames(c("estimate", "se", "p")) %>%
      setNames(c("(Intercept)_disp", "time_disp", "ontogenysmall_disp", "time:ontogenysmall_disp"))
  }
  
  TE <- getTimeEffect(fuu, transform = "response", ontogenylevel = "small")
  TE_col <- data.frame(timeeffect = c(TE$estimate, TE$SE, TE$p.value))
  
  OE <- getOntogenyEffect(fuu, transform = "response")
  OE_col <- data.frame(ontogenyeffect = c(OE$estimate, OE$SE, OE$p.value))
  
  isrelevant <- OE$signif & (TE$signif_rejneg | TE$signif_rejpos)
  direction <- if(isrelevant) { if(sign(TE$estimate) == sign(OE$estimate)) "e" else "o" } else ""

  # opposite <- sign(OE$estimate) != sign(TE$estimate) & OE$p.value < 0.05
  
  D <- getData(fuu)
  
  Years <- c(str_split_fixed(D$years,"–" ,n = 2), "")  # vector[3]!
  Years[1] <- paste0(Years[1], "–")
    
  
  # sd_ranef <- attr(s$varcor$cond$clusterid, "stddev")
  # n_clusters <- as.numeric(s$ngrps$cond)
  n_clusters_ES <- nrow(s$coefficients$cond) - 4
  
  
  R <- bind_cols(Titles, type = rownames(C_cond),
                 C_cond, C_disp,
                 OE_col, TE_col, direction = direction,
                 n = s$nobs,
                 n_clusters_ES = replace(n_clusters_ES, n_clusters_ES == 0, NA),
                 n_plots = D$n_plots,
                 yr_betweenAvg = D$yr_betweenAvg,
                 years = Years,
                 yr_max = D$yr_max,
                 year1 = D$year1,
                 year2 = D$year2
                 # n_clusters = ifelse(is.null(n_clusters), NA, n_clusters),
                 # sd_ranef = ifelse(is.null(sd_ranef), NA, sd_ranef)
  ) %>%
    mutate(type = fct_recode(type, estimate = 'Estimate',  se = 'Std. Error', p = "Pr(>|z|)"))
  
  rownames(R) <- NULL
  
  # time <- getEffectColumns(fuu, id)
  return(R)
}

Table <- lapply(fits, printModelRows) %>% bind_rows()

# # Set the first if arm to FALSE to enable the second. This way all estimates. also DE and ES will be tested.
# Check <- data.frame(plogis(Table$`(Intercept)_cond` + Table$ontogenysmall_cond) -
#                               plogis(Table$`(Intercept)_cond`), Table$ontogenyeffect, Table$inventory, Table$type) %>%
#   filter(Table.type == "estimate")
# # The error is on the response scale, 0.001 corresponds to 0.1 percent. This is probably numerics. Was worse for Spain with attime = 0 (insufficient centering)

# boxplot((Check[,1] - Check[,2]) ~ Check$Table.inventory)


saveRDS(Table, glue("Fits/{loadrunid}_Table_{modelname}_{loadrundate}.rds"))
