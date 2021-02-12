# The Library ---------------------------------------------------------
## Wrangling.
library(tidyverse) # ggplot2, purrr, tibble, dplyr, tidyr, stringr, readr, forcats
library(magrittr)
library(glue)
library(stringi)

## Mapping
library(sf)
library(raster)
library(eurostat)
library(elevatr)

## Plotting
library(ggplot2)
library(ggspatial)
library(cowplot)
library(gridGraphics)
library(gridExtra)
library(ggpubr)
library(ggfx)
library(extrafont)
library(colorspace)

## Tabulating.
library(stargazer)




# Set working directory ---------------------------------------------------
setwd(rprojroot::find_root(rprojroot::criteria$is_git_root))
inventoriespath <- '../Inventories/'


# Load theme and set color scales ----------------------------------------------------------
if (file.exists("Theme.R")) source("Theme.R") else theme_ockham <- function(rangeframemap) NULL
basecolor <- "#13223F"
inversebasecolor <- "#F7F8F9"
intermediatecolor <- "#A9AFB9"
watercolor <- "#33425F" # "#101620"
landcolor <- "#303030" # "#F5F6F7" # "#EFF0F2"
bordercolor <- "#5A5A5A" # "#F5F6F7" # "#EFF0F2"
backgroundcolor <- "#F5F6F7" # "#697285" "#A9AFB9" "#D4D7DC" "#EAEBEE" "#EFF0F2" "#F5F6F7" "#F7F8F9" are blends between basecolor and white


n_breaks_boxplot <- 4 

colorscale_cwb <- scale_color_gradient2(high = "#13225F", # "#082586", # bluer viridis violet
                                        mid = "white", ## "#FFDBFF", # similar to viridis green
                                        low = "#FF2010", #redder mennig ## "#FF8536", # lighter pumpkin
                                        midpoint = 0,
                                        aesthetics = c("colour", "fill"),
                                        breaks = scales::pretty_breaks(n = n_breaks_boxplot))


# scale_color_continuous(type = "viridis", direction = -1, aesthetics = c("colour", "fill"))

colorscale_gdd <- scale_color_gradient(low = "#13225F", # "#082586", # bluer viridis violet
                                       high = "#FF7036", # pumpkun-mennig
                                       aesthetics = c("colour", "fill"),
                                       breaks = scales::pretty_breaks(n = n_breaks_boxplot))

colorscale_direction_gdd <- scale_colour_discrete(type = c(basecolor, # green
                                                           "#FF4D26"  # mennig
), na.value = "transparent")

colorscale_direction <- scale_colour_discrete(type = c(intermediatecolor, basecolor), # intermediate bluegrey
                                              na.value = "transparent")

inventorycolor <-  c(ES = "#FF4D26", #mennig ## "#FF7526", # ES: pumpkin # "#FF5D26", # lighter small red
                     DE = "#5FEF8B", #mint ## "#16DB93", # DE: green # "#951147", # bordeaux in between
                     SE = basecolor) # "#382576") # SE: viridis violet

getColorscale_inv <- function(aesthetics, direction = 1, ...) scale_colour_discrete(type = unname(inventorycolor),
                                                                                    direction = direction,
                                                                                    aesthetics = aesthetics,
                                                                                    na.value = "transparent", ...) # special value for insignificant counts




# Load fit results -----------------------------------------------------------
loadrunid <- "DEF"
loadrundate <- "2020-11-02" # format(Sys.time(), '%Y-%m-%d')
modelname <- "glmmTMB_beta_ontogeny_time"

Table <- readRDS(glue("Fits/{loadrunid}_Table_{modelname}_{loadrundate}.rds"))  # Regression table for summary statistics

Effects <- readRDS(glue("Fits/{loadrunid}_Effects_{modelname}_{loadrundate}.rds")) %>%
  filter(env != "p_esdacc_01") %>%
  filter(env != "p_esdacc") %>%
  bind_cols(Table[match(.$id, Table$id), c("yr_betweenAvg", "years", "yr_max", "year1", "year2")]) %>% # crucial that match only extracts first occurrence!
  droplevels()



##########################################################################
# Wrangle plot values    ----------------------------------------------##
##########################################################################

isOutlier <- function(x) {
  x < quantile(x, 0.25) - (1 * IQR(x)) | x > quantile(x, 0.75) + (1 * IQR(x)) # usually (1.5 * IQR(x))
}

isExtreme <- function(x) {
  x == min(x) | x == max(x)
}

getAt <- function(targetcol, withincol, at) {
  i <- which(withincol == at)
  targetcol[i]
}

getM <- function(x, se) {
  m <- weighted.mean(x, 1/(se^2), na.rm = T)
  # m[is.nan(m)] <- 0
  return(m)
}

E <- Effects %>%
  # mutate(effect = fct_rev(effect)) %>% # for temporal to appear first
  mutate(enveff = interaction(env, effect), invenv = interaction(inventory, env, lex.order = T)) %>%
  mutate(inv = fct_recode(inventory, ES = "ES_IFN", DE = "DE_BWI", SE = "SE_NFI")) %>%
  mutate(tax = str_replace(tax, fixed("."), " "), tax_short = gsub("(?<=[A-Z])[a-z]*", ".", tax, perl = T)) %>%
  mutate(effectdirection = sign(estimate)) %>% # 1 or -1
  
  ## average estimates within inventory/env/effect
  group_by(inventory, env, effect) %>%
  mutate(m_inv.env.eff_pos = getM(ifelse(estimate > 0, estimate, NA), SE), m_inv.env.eff_neg = getM(ifelse(estimate < 0, estimate, NA), SE)) %>%
  mutate(displayname_enveff = ifelse(isExtreme(estimate), tax_short, NA)) %>%
  ungroup() %>%
  
  ## average estimates within env/effect/inventory
  group_by(env, effect, inventory) %>%
  mutate(m_env.eff.inv = getM(estimate, SE), m_env.eff.inv_pos = getM(ifelse(estimate > 0, estimate, NA), SE), m_env.eff.inv_neg = getM(ifelse(estimate < 0, estimate, NA), SE)) %>%
  mutate(displayname_inveff = ifelse(isExtreme(estimate), tax_short, NA)) %>%
  ungroup() %>%
  
  ## average estimates within tax/env/effect
  group_by(tax, env, effect) %>%
  mutate(m_tax.env.eff = getM(estimate, SE), m_tax.env.eff_pos = getM(ifelse(estimate > 0, estimate, NA), SE), m_tax.env.eff_neg = getM(ifelse(estimate < 0, estimate, NA), SE)) %>%
  ungroup() %>%
  
  ##
  group_by(tax, env) %>%
  mutate(m_tax.env.ont = getM(ifelse(effect == "ontogeny", estimate, NA), SE)) %>%
  mutate(m_abs_signif_tax.env.ont = getM(ifelse(effect == "ontogeny" & signif, abs(estimate), NA), SE), m_abs_signif_tax.env.time = getM(ifelse(effect == "time" & signif, abs(estimate), NA), SE)) %>%
  ## order by 
  mutate(order_ont = if_else(m_abs_signif_tax.env.ont > 0, m_abs_signif_tax.env.ont + 1000, getM(ifelse(effect == "ontogeny", abs(estimate), NA), SE), missing = 0),
         order_time = if_else(m_abs_signif_tax.env.time > 0, m_abs_signif_tax.env.time + 1000, getM(ifelse(effect == "time", abs(estimate), NA), SE), missing = 0)) %>%
  ungroup() %>%
  
  ## effects significantly different and differently directed?
  ## getAt gets the position condition-matched within one vector from another within the group
  group_by(inv, tax, env) %>%
  ## Are the direction variables even used?
  mutate(ontdirection = getAt(targetcol = effectdirection, withincol = effect, at = "ontogeny")) %>%
  mutate(diffdirection = signif_diff & !(var(effectdirection, na.rm = T) == 0)) %>%
  mutate(diffdirection = ifelse(diffdirection, "!", NA)) %>%   ## TRUE become !, everything else NA
  
  ## See below: NA inflation for plotting reasons (NA -> invisible)
  mutate(tempsignif.onesided = getAt(signif_rejpos, effect, "time") | getAt(signif_rejneg, effect, "time")) %>%
  mutate(ontsignif = getAt(signif, effect, "ontogeny") | getAt(signif, effect, "ontogeny")) %>%
  mutate(ontsignif.diffdir = (effect == "ontogeny") & signif & !(var(effectdirection, na.rm = T) == 0)) %>%  ## significant and differently directed?
  mutate(ontsignif.eqdir = (effect == "time") & ontsignif & (var(effectdirection, na.rm = T) == 0)) %>%  ## significant and equal directed? (effect == "time") & 

  mutate(ontsignif.eqdir.tempsignif = ontsignif.eqdir & tempsignif.onesided & effect == "time") %>% ## This puts it into the time rows!
  mutate(ontsignif.diffdir.tempsignif = ontsignif.diffdir & tempsignif.onesided & effect == "ontogeny") %>% ## This puts it into the ontogeny rows!
  
  mutate(ontsignif.eqdir.tempsignif.plotleft = getAt(ontsignif.eqdir, effect, "time") & getAt(tempsignif.onesided, effect, "time") & effect == "ontogeny") %>% ## This puts it into the ontogeny rows for plotting
  
  
  mutate_at(vars(starts_with(c("tempsignif.", "ontsignif."))), function(x) replace(x, x == FALSE, NA)) %>%
  ## No NA inflaction here!
  mutate(ontsignif.tempsignif.onesided = ontsignif & getAt(tempsignif.onesided, withincol = effect, "ontogeny")) %>%
  
  ungroup() %>%
  
  ## only plot years on the time effect side
  mutate(yr_betweenAvg.plot = ifelse(effect == "time", yr_betweenAvg, NA)) %>%
  
  ## order manually south to north (for north to south plotting)
  mutate(inventory = factor(inventory, levels = c("ES_IFN", "DE_BWI", "SE_NFI")),
         inv = factor(inv, levels = c("ES", "DE", "SE")),
         invenv = factor(invenv, levels = c("ES_IFN.cwbYear_aclim_01", "ES_IFN.gdd0_aclim_01", 
                                            "DE_BWI.cwbYear_aclim_01", "DE_BWI.gdd0_aclim_01", 
                                            "SE_NFI.cwbYear_aclim_01", "SE_NFI.gdd0_aclim_01")))

E_gdd <- filter(E, env == "gdd0_aclim_01") %>% droplevels()
E_cwb <- filter(E, env == "cwbYear_aclim_01") %>% droplevels()

E_count <- E %>%
  add_count(env, effect, effectdirection, wt = signif*effectdirection,  name = "n_signif") %>%
  add_count(env, effect, effectdirection, wt = effectdirection, name = "n_total") %>%
  mutate(n_add = n_total - n_signif) %>%
  mutate(inv = fct_expand(inv, "ns")) %>% # This adds a fourth inventory named 'ns'
  complete(nesting(env, effect, effectdirection), inv) %>% # this also completes one level in SE
  filter(is.na(n_total) & inv != "ns") %>%
  group_by(env, effect, effectdirection, n_total, n_signif, n_add, inv) %>%
  tally(wt = signif*effectdirection,  name = "n_signif_inventory") %>%
  ungroup()

# hacky base
E_count[is.na(E_count$n_add), "n_signif_inventory"] <- E_count[E_count$inv == "SE", "n_add"]
E_count[is.na(E_count$n_add), "n_total"] <- E_count[E_count$inv == "SE", "n_total"]
E_count[E_count$inv == "ns", "inv"] <- NA


E_count_gdd <- filter(E_count, env == "gdd0_aclim_01")
E_count_cwb <- filter(E_count, env == "cwbYear_aclim_01")



#######################################################################
# Map ------------------------------------------------------------------##
##########################################################################

# Load plot points --------------------------------------------------------
equalearthproj <- "+proj=laea +x_0=0 +y_0=0 +lon_0=0 +lat_0=0" # Lambert

plotid_pres <- readRDS("Data/taxtables_pres_thresholdsubset.rds") %>%
  bind_rows() %$%
  plotid %>%
  unique()
Plots <- readRDS(glue("{inventoriespath}/Combined data/Env.rds")) %>%
  filter(plotid %in% plotid_pres) %>%
  unique() %>%
  st_transform(equalearthproj) %>%
  mutate(country = fct_recode(inventory, SE = "SE_NFI", DE = "DE_BWI", ES = "ES_IFN")) %>%
  mutate(country = fct_relevel(country, levels = c("SE", "DE", "ES")))

Envrange <- data.frame(gdd = with(Plots, c(min(gdd0_aclim), max(gdd0_aclim))),
                       cwb = with(Plots, c(min(cwbYear_aclim), max(cwbYear_aclim))))
print(Envrange)

adminid <- c("SE", "DE", "ES")

# Areas <- eurostat::get_eurostat_geospatial(resolution = "1", nuts_level = "0") %>%
#   bind_rows(eurostat::get_eurostat_geospatial(resolution = "1", nuts_level = "1")) %>%
#   dplyr::filter(id %in% !!adminid) %>%
#   st_transform(equalearthproj) %>%
#   sf::st_buffer(dist = 0) %>% # !
#   st_crop(xmin = -1200000,
#           xmax = 1680000,
#           ymin = 3800000, # Crop southern islands.
#           ymax = 7500000)
# saveRDS(Areas, "Data/Plotareas_SE_DE_ES_sf.rds")

# Background <- eurostat::get_eurostat_geospatial(resolution = "1", nuts_level = "0") %>%
#   st_transform(equalearthproj) %>%
#   dplyr::filter(!(id %in% !!c("IS"))) %>%
#   sf::st_buffer(dist = 0) %>% # fixes intersection problem
#   st_crop(xmin = -1200000,
#           xmax = 1680000,
#           ymin = 3800000, # Crop southern islands.
#           ymax = 7500000)
# saveRDS(Background, "Data/Plotareas_background_sf.rds")

Areas <- readRDS("Data/Plotareas_SE_DE_ES_sf.rds")
Background <- readRDS("Data/Plotareas_background_sf.rds")

# plot(st_geometry(Areas), border = "grey")
# plot(Plots["gdd0_aclim"], cex = 0.01, pal = colorRampPalette(c("#FF4D26", nicecolor["lightblue"])), add = T)


# Map ---------------------------------------------------------------------

mapPlotsLight <- function(env = "gdd0_aclim",
                     size = 1e-100, # there seems to be some lower threshold here, under which no size can be provided.
                     title = NULL,
                     colorscale = scale_color_continuous(type = "viridis", direction = -1)) {
  map <- ggplot(Plots) +
    theme(text = element_text(color = basecolor, family = ".SF Compact Text", size = 14)) +
    geom_sf(data = Background, fill = backgroundcolor, colour = backgroundcolor) +
    geom_sf(data = Areas, colour = basecolor, fill = "transparent") +
    geom_sf(aes_string(color = env), size = size) +
    guides(color = guide_colourbar(barwidth = 0.4, barheight = 8, draw.ulim = TRUE, draw.llim = TRUE, frame.linewidth = 0, nbin = 100)) +
    theme(legend.position = c(0.12, 0.8)) +
    { if (!is.null(title)) ggtitle(title)} +
    colorscale +
    theme_empty()
  return(map)
}


mapPlots <- function(env = "gdd0_aclim",
                     size = 1e-100, # there seems to be some lower threshold here, under which no size can be provided.
                     title = NULL,
                     legendtitle = NULL,
                     mapaccidentiae = T,
                     plotcountrylabels = T,
                     colorscale = scale_color_continuous(type = "viridis", direction = -1)
                     ) {
  
  
  ## Map
  map <- ggplot(Plots) +
    # theme(plot.background = element_rect(fill = watercolor)) +
    theme(legend.background = element_rect(fill = NA)) +
    theme(text = element_text(color = basecolor, family = ".SF Compact Text", size = 14), plot.title = element_text(hjust = 0.5)) +
    
    with_blur(geom_sf(data = Background, fill = watercolor, colour = watercolor)[[1]], # index extracts the first item, a layer object
              sigma = 9) + 
    geom_sf(data = Background, fill = landcolor, colour = landcolor) +
    geom_sf(data = Areas, colour = bordercolor, fill = "transparent", lwd = 0.9) +
    { if (plotcountrylabels) annotate("text", x = c(-5, 5.9, 18.2), y = c(44.7, 54.8, 56), label = c("ES", "DE", "SE"), col = inventorycolor, size = 5.3) } +
    theme(axis.title = element_blank()) + # Somehow this is necessary again for annotate.
    
    geom_sf(aes_string(color = env), size = size) +
    
    guides(color = guide_colourbar(barwidth = 0.45, barheight = 10, ticks.linewidth = 0.8, draw.ulim = TRUE, draw.llim = TRUE, frame.linewidth = 1, nbin = 100, title = legendtitle)) +
    theme(legend.position = c(0.1, 0.8)) +
    { if (!is.null(title)) ggtitle(title)} +
    colorscale +
    theme_empty() +
    
    { if(mapaccidentiae) annotation_scale(location = "br", width_hint = 0.5, style = "bar",
                     bar_cols = c(landcolor, "white"), line_col = landcolor, text_col = basecolor, text_cex = 1.15, text_family = ".SF Compact Text")} +
    { if(mapaccidentiae) annotation_north_arrow(location = "bl", which_north = "true",
                           pad_x = unit(0.2, "cm"),
                           style = north_arrow_minimal(line_col = landcolor, text_col = basecolor, fill = c(landcolor), text_size = 17, text_family = ".SF Compact Text")) }
  
  ## Boxplot
  countryboxplot <- ggplot(Plots, aes_string(x = "country", y = env)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = n_breaks_boxplot)) +
    geom_boxplot(fill = NA, outlier.size = 0.5, outlier.shape = NA) +
    theme_ockham()
  
  map <- ggdraw(map) +
    draw_plot(countryboxplot, x = 0.15, y = 0.55, width = 0.4, height = 0.35)
  
  return(map)
}


# map_gdd <- mapPlotsLight("gdd0_aclim", colorscale = colorscale_gdd)
# map_cwb <- mapPlotsLight("cwbYear_aclim", colorscale = colorscale_cwb)

map_gdd_big <- mapPlots("gdd0_aclim", title = "Heat sum (growing degree days)", colorscale = colorscale_gdd, mapaccidentiae = F, legendtitle = "heat sum [°C d]")
map_cwb_big <- mapPlots("cwbYear_aclim", title = "Water availability (climatic water balance)", colorscale = colorscale_cwb, plotcountrylabels = F, legendtitle = "water availability [mm]")


mapgrid_gdd <- cowplot::plot_grid(map_gdd_big, map_cwb_big, ncol = 2, labels = "AUTO", label_size = 16)
cowplot::save_plot(glue("Publishing/Plots/Map_gdd_cwb.pdf"), mapgrid_gdd, base_height = 9, base_asp = 1.4, device = cairo_pdf) # cairo_pdf will embed fonts!


######################################################################################
# Effect comparision ----------------------------------------------------------------##
######################################################################################


# effectplot_allinone <- ggplot(E,
#                      mapping = aes(x = enveff, y = estimate, ymin = lower.CL, ymax = upper.CL, group = effect, col = env)) +
#   theme_ockham(xlabangle = 60, rangeframemap = aes(y = NULL)) + ## the theme.
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
#   facet_wrap(~ inventory) + ## facet layout
# 
#   ## effect estimates
#   geom_pointrange(position = position_jitterdodge(jitter.width = 0.01, dodge.width = c(0.1), seed = 1), alpha = 0.15) +
#   geom_hline(yintercept=0, linetype="dotted") +
#   
#   ## weighted average effects 
#   geom_segment(aes(x = enveff, xend = enveff, y = 0, yend = m_inv.env.eff_pos),
#                arrow = arrow(angle = 20, ends = "last", length = unit(6, "pt")), size = 0.5, position = position_nudge(0.08)) +
#   geom_segment(aes(x = enveff, xend = enveff, y = 0, yend = m_inv.env.eff_neg),
#                arrow = arrow(angle = 20, ends = "last", length = unit(6, "pt")), size = 0.5, position = position_nudge(-0.08)) +
#   ## additional labeling
#   geom_text(aes(label = displayname_enveff), na.rm = TRUE, nudge_y = 0.002, cex = 2)
#   # scale_discrete_manual(values = c("#FF4D26", "#00C3FF"), c("colour", "fill"))
# 
# effectplot_allinone


# effectplot_ontogeny <- ggplot(filter(E, effect == "ontogeny"), mapping = aes(x = invenv, y = estimate, ymin = lower.CL, ymax = upper.CL, col = env)) +
#   theme_ockham(rangeframemap = aes(x = NULL)) + ## the theme.
#   scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(-0.09, 0.09)) +
#   
#   ## effect estimates
#   geom_pointrange(position = position_jitterdodge(jitter.width = 1, dodge.width = c(0.001), seed = 1), alpha = 0.15) +
#   geom_hline(yintercept=0, linetype="dotted") +
#   
#   ## weighted average effects 
#   geom_segment(aes(x = invenv, xend = invenv, y = 0, yend = m_inv.env.eff_pos),
#                arrow = arrow(angle = 20, ends = "last", length = unit(6, "pt")), size = 0.5, position = position_nudge(0.05)) +
#   geom_segment(aes(x = invenv, xend = invenv, y = 0, yend = m_inv.env.eff_neg),
#                arrow = arrow(angle = 20, ends = "last", length = unit(6, "pt")), size = 0.5, position = position_nudge(-0.05)) +
#   ## additional labeling
#   geom_text(aes(label = displayname_enveff), na.rm = TRUE, nudge_x = 0.02, cex = 2.5) +
#   # scale_discrete_manual(values = c("#FF4D26", "#00C3FF"), c("colour", "fill")) +
#   coord_flip() +
#   ggtitle("ontogeny effect")



plotAverage <- function(E_sub, rangeframemap, ...) {
  
  effectplot <- ggplot(E_sub, mapping = aes(x = inv, y = estimate, ymin = lower.CL, ymax = upper.CL, col = estimate > 0)) +
    
    ## setting up the background grid
    geom_tile(aes(height = Inf, fill = inv), lwd = 0, colour = "transparent") +
    scale_fill_manual(values = rep_len(c(backgroundcolor, "transparent"), nrow(E_sub)), guide = "none") + # "#697285" "#A9AFB9" "#D4D7DC" "#EAEBEE" "#EFF0F2" "#F5F6F7" "#F7F8F9" are blends between basecolor and white
    geom_hline(yintercept = 0, linetype= "dotted") +
    
    ## effect estimates
    # geom_pointrange(position = position_jitter(width = 0.2, seed = 10), alpha = 0.2, size = 1, fatten = 0.9) + # fatten is multiplicative in/decreases point compared to bar
    geom_point(size = 1.5, alpha = 0.2, shape = 16, position = position_jitter(width = 0.2, seed = 10)) +
    geom_linerange(size = 0.5, alpha = 0.2, position = position_jitter(width = 0.2, seed = 10)) +
    
    ## weighted average effects 
    geom_segment(aes(x = inv, xend = inv, y = 0, yend = m_env.eff.inv_pos, col = T),
                 arrow = arrow(angle = 20, ends = "last", length = unit(8, "pt")), size = 0.8, position = position_nudge(0.07)) +
    geom_segment(aes(x = inv, xend = inv, y = 0, yend = m_env.eff.inv_neg, col = F),
                 arrow = arrow(angle = 20, ends = "last", length = unit(8, "pt")), size = 0.8, position = position_nudge(-0.07)) +
    colorscale_direction_gdd + # color scale for arrows etc.
    
    ## additional labeling
    # geom_text(aes(label = displayname_inveff), na.rm = TRUE, nudge_x = 0.02, cex = 2.5) +
    
    ## axes and layout
    xlab(NULL) + ylab(NULL) +
    ggtitle("Title") +
    coord_flip() +
    facet_wrap(~ effect, labeller = as_labeller(function(f) paste(ifelse(f == "ontogeny", "life stage", f), "shift"))) + ## facet layout
    
    theme_ockham(rangeframemap = aes(x = NULL),
                 panel.spacing.x = unit(-8, "lines"),
                 axis.ticks.y = element_blank(),
                 ...
    ) + ## the theme.
    
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) # limits = c(-0.09, 0.09)
  
  grid <- cowplot::plot_grid(rangeframemap, effectplot, nrow = 1, rel_widths = c(3, 3.3))
  return(grid)
}


averagegrid_gdd <- plotAverage(E_gdd, map_gdd, legend = "bottom")
cowplot::save_plot(glue("Publishing/Plots/{loadrunid}_Average_gdd.pdf"), averagegrid_gdd, base_height = 7, base_asp = 2, device = cairo_pdf, bg = "transparent") # cairo_pdf will embed fonts!

# averagegrid_cwb <- plotAverage(E_cwb, map_cwb, legend = "bottom")
# cowplot::save_plot(glue("Publishing/Plots/{loadrunid}_Average_cwb.pdf"), averagegrid_cwb, base_height = 7, base_asp = 2, device = cairo_pdf) # cairo_pdf will embed fonts!


######################################################################################
# Species comparision --------------------------------------------------------------##
######################################################################################


plotTaxa <- function(E_sub, E_count, orderby = quo(order_ont),
                     plotperiods = F, plotaveragebydirection = F,
                     ...) { # ... passed on to theme_ockham() which passes on to theme
  E_sub <- mutate(E_sub, diffdirection = ifelse(effect == "time", NA, diffdirection)) %>% # empty the "!" labels in the time facet
    mutate(tax = fct_reorder(tax, !!orderby))
  
  dodgymcdodgeface <- position_dodge2(width = 0.7, preserve = 'single')
  panelspacing <- unit(-1, "lines")
  
  upper.CL.isveryhigh <- E_sub$upper.CL > 0.065
  E_sub$upper.CL[upper.CL.isveryhigh] <- 0.065
  
  addinfoypos <- max(E_sub$upper.CL) + 0.005
  
  
  taxplot <- ggplot(E_sub, 
                    aes(tax, # Ordered by the weighted mean of the absolute ontogeny effect size, given significance.
                        estimate,
                        color = inv,
                        group = interaction(tax, inv),
                        ymin = lower.CL,
                        ymax = upper.CL)
                    ) +
    
    #### This is just the background grid
    geom_tile(aes(height = Inf, fill = tax), lwd = 0, colour = "transparent") +
    scale_fill_manual(values = rep_len(c(backgroundcolor, "transparent"), nrow(E_sub)), guide = "none") +
    geom_hline(yintercept=c(0), linetype = "dotted") +
    geom_hline(yintercept=c(-0.05, 0.05), linetype = "dotted", color = "#D4D7DC") +
    
    #### These are the point/lineranges
    # geom_pointrange(aes(stroke = signif, alpha = signif), position = dodgymcdodgeface, fatten = 3) +
    # scale_shape_tremmel() +
    geom_point(aes(shape = signif), position = dodgymcdodgeface, size = 2.2, stroke = 0.5, alpha = 1) + # aes(size = signif, alpha = signif)
    # scale_size_manual(values = c(2.5, 2)) + # larger for nonsignificant
    geom_linerange(aes(alpha = signif), position = dodgymcdodgeface) +
    scale_shape_manual(values = c(1, 16)) + # Solid for significant
    scale_alpha_discrete(range = c(0.42, 1)) + 
    # scale_size_discrete(range = c(1, 2)) +
    getColorscale_inv(aesthetics = c("colour"), guide = "none") +
    
    #### Mark different/equal direction
    # geom_text(aes(y = max(upper.CL) - (as.numeric(inv)-1)*0.004, x = tax, label = diffdirection), size = 3, alpha = 1) + # significantly different and different direction
    # geom_point(aes(y = max(upper.CL) + ontsignif.diffdir - 1.008, x = tax, alpha = ontsignif.diffdir), size = 1.5, shape = 8, position = dodgymcdodgeface, guide = "none") + # significantly different and different direction
    
    # geom_point(aes(y = max(upper.CL) + ontsignif.diffdir.tempsignif - 1.015, x = tax, alpha = ontsignif.diffdir.tempsignif), size = 2.5, shape = 8, position = dodgymcdodgeface, guide = "none") + # different direction
    # geom_point(aes(y = max(upper.CL) + ontsignif.eqdir.tempsignif.plotleft - 1.005, x = tax, alpha = ontsignif.eqdir.tempsignif.plotleft), size = 2.5, shape = 8, position = dodgymcdodgeface, guide = "none") + # equal direction
  
    geom_text(aes(y = addinfoypos + 0.008 + ontsignif.diffdir.tempsignif - 1, x = tax, alpha = ontsignif.diffdir.tempsignif), size = 3, label = "o", position = dodgymcdodgeface, guide = "none") + # different direction
    geom_text(aes(y = addinfoypos + 0.016 + ontsignif.eqdir.tempsignif.plotleft - 1, x = tax, alpha = ontsignif.eqdir.tempsignif.plotleft), size = 3, label = "e", position = dodgymcdodgeface, guide = "none") + # equal direction
    # geom_hline(yintercept=c(addinfoypos + 0.03), size = 0.1) +
    
    
    #### Add time Bars
    { if(plotperiods)
    geom_hline(yintercept = c(addinfoypos), size = 0.1) +
    geom_hline(yintercept = addinfoypos + 0.0008*c(10, 20), linetype = "dotted", color = "#D4D7DC") + # dotted lines for 10, 20 years
    geom_linerange(aes(ymin = addinfoypos, ymax = addinfoypos + yr_betweenAvg.plot*0.0008, x = tax), size = 1.1, alpha = 0.6, position = dodgymcdodgeface, guide = "none") # different direction
    }  +
    
    
    #### Facets
    facet_wrap(~ effect, labeller = as_labeller(function(f) fct_recode(f, "juvenile divergence" = "ontogeny", "temporal shift of juveniles" = "time"))) +
    
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    
    xlab(NULL) + ylab(NULL) +
    # ggtitle("Title") +
    theme_ockham(rangeframemap = aes(x = NULL), # rangeframemap is the aes for the range axes
                 axis.ticks.y = element_blank(),
                 axis.text.y = element_text(face = "italic"),
                 # legend = 'none',
                 panel.spacing.x = panelspacing, ...) + # -1 , ...
    coord_flip()
  
  
  yplotrange <- ggplot_build(taxplot)$layout$panel_scales_y[[1]]$range$range
  
  ## old version with stacked bars
  # taxcountplot <- ggplot(E_count, aes(x = effectdirection, y = n_signif_inventory, size = abs(n_total), fill = inv)) +
  #     theme_ockham(axes = "none") +
  #     facet_wrap(~ effect, labeller = as_labeller(function(f) glue("count of significant {f} effects by direction
  #                                                                THIS WILL BE ALIGNED BY 0 WITH THE ABOVE SCALE, LABELS REPOSITIONED"))) +
  #     geom_col(position = position_stack(reverse = T), lwd = 0.28, col = scales::alpha(basecolor, 0.6)) + # col is only for frame
  #     getColorscale_inv(aesthetics = c("fill")) +
  # 
  #     geom_text(aes(label = glue("{n_signif_inventory}
  #                                {inv}")), position = position_stack(vjust = .5, reverse = T), size = 2.5) +
  #     geom_text(aes(x = effectdirection, y = n_total + 4*effectdirection, label = glue("total
  #                                                                                    {n_total}")), size = 2.5) +
  #     coord_flip()
  #   
  
  
  E_count <- E_sub %>%
    ## 1, separately make all counts that are conditional on ontsignif (plot row D)
    filter(ontsignif) %>% # filter only groups (models including the time effect), where there is a significant ontigeny effects
    mutate(signif = ontsignif.eqdir | ontsignif.diffdir) %>% # both cases get a count ...
    mutate(effect = case_when(ontsignif.eqdir ~ "ontogeny_eq",
                              ontsignif.diffdir ~ "ontogeny_diff")) %>% # ... but get a different effect name (== facet)
    ## Classifying into three types + NA
    mutate(ontsignif.type = case_when(ontsignif.eqdir.tempsignif ~ "ontsignif.eqdir.tempsignif",
                                      ontsignif.diffdir.tempsignif ~ "ontsignif.diffdir.tempsignif",
                                      ontsignif ~ "ontsignif")) %>% # 3 groups: whether ontsignif is either also a, b, or not
    

    ## 2 bind all other counts
    bind_rows(E_sub) %>%
    mutate(inv = replace(as.character(inv), !signif | is.na(signif), NA)) %>%
    filter(!is.na(inv)) %>%
    
    mutate(tempsignif.star = effect %in% c("ontogeny_diff", "ontogeny_eq") & ontsignif.tempsignif.onesided) %>% ## includes NAs for plotting
    mutate(tempsignif.star.tf = replace(tempsignif.star, is.na(tempsignif.star), FALSE)) %>% # has only TRUE and FALSE for group plotting
    
    ## ! to add a fourth type classify
    mutate(ontsignif.type = case_when(is.na(ontsignif.type) & signif ~ "signif",
                                      # for the types that had been there before
                                      TRUE ~ ontsignif.type)) %>%
    # mutate(ontsignif.type = factor(ontsignif.type, levels = c("signif", "ontsignif", "ontsignif.diffdir.tempsignif", "ontsignif.eqdir.tempsignif"))) %>%
    

    
    add_count(inv, effect, effectdirection, name = "count") %>%
    mutate(inv = factor(inv, levels = c("ES", "DE", "SE")), effect = factor(effect, levels = c("ontogeny", "time", "ontogeny_diff", "ontogeny_eq"))) %>%
    # alternatively for temporal first: # mutate(inv = factor(inv, levels = c("ES", "DE", "SE")), effect = factor(effect, levels = c("time", "ontogeny", "ontogeny_eq", "ontogeny_diff"))) %>%
    arrange(ontsignif.type) # order the data.frame, because this order is how position_stack stacks the points
  
  # mutate(dotplotsymbol = fct_recode(effect, "one" = "ontogeny", "one" = "time", "eight" = "ontogeny_diff"))
  #   
  # 
  # ## another version with stacked bars
  # taxcountplot <- ggplot(E_count, aes(x = inv, y = effectdirection, color = inv, fill = inv), group = inv) +
  #     theme_ockham(axes = "none") +
  #     theme(strip.text.x = element_blank()) +
  #     facet_wrap(~ effect, nrow = 1) + # labeller = as_labeller(function(f) glue("count of significant {f} effects by direction, THIS WILL BE ALIGNED BY 0 WITH THE ABOVE SCALE, LABELS REPOSITIONED"))
  #     geom_dotplot(stackgroups = F, binpositions = "bygroup", stackratio = 1.4, binwidth = 0.2, method = "histodot", drop = T) +
  #     getColorscale_inv(aesthetics = c("colour", "fill"), direction = 1) +
  # 
  #     # geom_text(stat = 'count', aes(label = ..count.., y = 0.5), hjust = -7) +
  #     geom_text(aes(y = -3 , label = inv)) +
  #   
  #     # geom_text(aes(x = effectdirection, y = n_total + 4*effectdirection, label = glue("total
  #     #                                                                                {n_total}")), size = 2.5) +
  # 
  #    coord_flip()
  
  
  taxcountplot <- ggplot(E_count, 
                         aes(inv,
                             estimate, # this is just here for scale adjustment
                             color = inv,
                             shape = ontsignif.type, # tempsignif.star.tf, # effect # ontsignif.type
                             label = inv,
                             group = tax)) +
    
    ## This is just the background grid and everything else from the top plot
    facet_wrap(~ effect, nrow = 2, labeller = as_labeller(function(f) fct_recode(f, "no. of significant juvenile divergences" = "ontogeny",
                                                                                 "no. of significant temporal shifts" = "time",
                                                                                 "... with opposite temporal shift" = "ontogeny_diff",
                                                                                 "... with equal-direction temporal shift" = "ontogeny_eq"))) +
    geom_hline(yintercept=c(0), linetype = "dotted") +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = yplotrange) +
    xlab(NULL) + ylab(NULL) +
    theme_ockham(axes = "none",
                 rangeframemap = aes(x = NULL, y = NULL), # rangeframemap is the aes for the range axes
                 axis.text.y = element_text(colour = inventorycolor), # recode(as.character(E_count$inv), SE = "#482576", DE = "#1E9C89", ES = "#FDE725"))
                 panel.spacing.x = panelspacing) + # -1
    coord_flip() +
    
    geom_point(mapping = aes(x = inv,
                             y = 0.006 * effectdirection,
                             ), position = position_stack(vjust = 0.5), size = 2.4, stroke = 0.7) +
    scale_shape_manual(values = c(1, 16, 16, 16)) + # hollow for ontsignif, solid for signif, ontsignif.diffdir.tempsignif, and ...eqdir
    # this is a shape scale for three levels im ontsignif.type: # scale_shape_manual(values = c(16, 4, 3)) + # look at shape in aes
    geom_text(aes(y = count * 0.5 * 0.006 * effectdirection, label = count), position = position_nudge(x = 0.5), size = 3) +
    getColorscale_inv(aesthetics = c("colour", "fill"))
  
  
  if (plotaveragebydirection) {
    taxavgplot <- ggplot(E_sub, 
                         aes(inv,
                             estimate, # this is just here for scale adjustment
                             color = inv,
                             shape = effect,
                             label = inv,
                             group = inv)) +
      
      ## This is just the background grid and everything else from the top plot
      facet_wrap(~ effect, nrow = 1, labeller = as_labeller(function(f) fct_recode(f, "average juvenile divergence" = "ontogeny", "average temporal shift" = "time"))) +
      geom_hline(yintercept=c(0), linetype = "dotted") +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = yplotrange) + # c(min(E_sub$lower.CL), max(E_sub$upper.CL))
      xlab(NULL) + ylab(NULL) +
      theme_ockham(axes = "range",
                   axis.ticks.y = element_blank(),
                   rangeframemap = aes(x = m_env.eff.inv, y = m_env.eff.inv), # map is the aes for the range axes
                   axis.text.y = element_text(colour = inventorycolor), # recode(as.character(E_count$inv), SE = "#482576", DE = "#1E9C89", ES = "#FDE725"))
                   panel.spacing.x = panelspacing) + # -1
      coord_flip() +
      
      ## weighted average effects 
      geom_segment(aes(x = inv, xend = inv, y = 0, yend = m_env.eff.inv_pos, col = inv, alpha = T), # positive drection
                   arrow = arrow(angle = 20, ends = "last", length = unit(6, "pt")), size = 0.7, position = position_nudge(0.1)) +
      geom_segment(aes(x = inv, xend = inv, y = 0, yend = m_env.eff.inv_neg, col = inv, alpha = F), # negative drection
                   arrow = arrow(angle = 20, ends = "last", length = unit(6, "pt")), size = 0.7, position = position_nudge(-0.1)) +
      
      # geom_text(aes(y = count * 0.5 * 0.006 * effectdirection, label = count), position = position_nudge(x = 0.5), size = 3) +
      getColorscale_inv(aesthetics = "color")
    
  } else {
    
    ## uses m_env.eff.inv, the average of all effects, weighted by 1/SE
    taxavgplot <- ggplot(E_sub, 
                         aes(inv, # Ordered by the weighted mean of the absolute ontogeny effect size, given significance.
                             m_env.eff.inv, # this is just here for scale adjustment
                             color = inv,
                             shape = effect,
                             label = inv,
                             group = inv)) +
      
      ## This is just the background grid and everything else from the top plot
      facet_wrap(~ effect, nrow = 1, labeller = as_labeller(function(f) fct_recode(f, "average juvenile divergence" = "ontogeny", "average temporal shift" = "time"))) +
      geom_hline(yintercept=c(0), linetype = "dotted") +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(-0.015, 0.015)) + # c(min(E_sub$lower.CL), max(E_sub$upper.CL))
      xlab(NULL) + ylab(NULL) +
      theme_ockham(axes = "range",
                   axis.ticks.y = element_blank(),
                   axis.text.y = element_text(colour = inventorycolor), # recode(as.character(E_count$inv), SE = "#482576", DE = "#1E9C89", ES = "#FDE725"))
                   panel.spacing.x = panelspacing) + # -1
      coord_flip() +
      
      ## weighted average effects 
      geom_segment(aes(x = inv, xend = inv, y = 0, yend = m_env.eff.inv, col = inv), # positive drection
                   arrow = arrow(angle = 20, ends = "last", length = unit(6, "pt")), size = 0.7, position = position_nudge(0.1)) +
      geom_text(aes(y = m_env.eff.inv, label = percent(m_env.eff.inv, accuracy = 0.1)), position = position_nudge(x = 0.5), size = 3) +
      
      
      # geom_text(aes(y = count * 0.5 * 0.006 * effectdirection, label = count), position = position_nudge(x = 0.5), size = 3) +
      getColorscale_inv(aesthetics = "color")
  }

  
  grid <- cowplot::plot_grid(taxplot, taxavgplot, taxcountplot, nrow = 3, align = "v", axis = "blrt", rel_heights = c(9.1, 1.37, 2.5),
                             labels = "AUTO", label_size = 14)
  
  return(grid)
}

## CAUTION! Letters indicating direction might be covered by facets layered on top. Keep in mind during post-editing.

taxgrid_gdd <- plotTaxa(E_gdd, E_count_gdd) # geom_text warnings are ok
cowplot::save_plot(glue("Publishing/Plots/{loadrunid}_Taxa_gdd.pdf"), taxgrid_gdd, base_height = 11.5, base_asp = 0.8, device = cairo_pdf, bg = "transparent") # cairo_pdf will embed fonts!
# embed_fonts(glue("Publishing/Plots/{loadrunid}_Taxa_gdd.pdf"))

taxgrid_cwb <- plotTaxa(E_cwb, E_count_cwb) # geom_text warnings are ok
cowplot::save_plot(glue("Publishing/Plots/{loadrunid}_Taxa_cwb.pdf"), taxgrid_cwb, base_height = 11.5, base_asp = 0.8, device = cairo_pdf, bg = "transparent")

taxgrid_gdd_bytime <- plotTaxa(E_gdd, E_count_gdd, orderby = quo(order_time)) # geom_text warnings are ok
cowplot::save_plot(glue("Publishing/Plots/{loadrunid}_Taxa_gdd_bytime.pdf"), taxgrid_gdd_bytime, base_height = 11.5, base_asp = 0.8, device = cairo_pdf, bg = "transparent") # cairo_pdf will embed fonts!
# embed_fonts(glue("Publishing/Plots/{loadrunid}_Taxa_gdd.pdf"))

taxgrid_cwb_bytime <- plotTaxa(E_cwb, E_count_cwb, orderby = quo(order_time)) # geom_text warnings are ok
cowplot::save_plot(glue("Publishing/Plots/{loadrunid}_Taxa_cwb_bytime.pdf"), taxgrid_cwb_bytime, base_height = 11.5, base_asp = 0.8, device = cairo_pdf, bg = "transparent")


######################################################################################
# Regression table -----------------------------------------------------------------##
######################################################################################
library(kableExtra)

keepFirstRow <- function(string) if_else(row_number() == 1, as.character(string), '')

pasteFormatted <- function(num, t) {
  num_char <- num
  num_char[t == "p"] <- symnum(num[t == "p"], corr = FALSE, na = FALSE, cutpoints = c(0, 0.05, 1), symbols = c("*", " "))
  num_char[t == "se"] <- glue("({sprintf('%.3f', num[t == 'se'])})")
  num_char[t == "estimate"] <- sprintf('%.4f', num[t == "estimate"])
  num_char[t == "percent"] <- sprintf('%.2f%s', num[t == "percent"] * 100, "%")
  return(num_char)
}

Table <- readRDS(glue("Fits/{loadrunid}_Table_{modelname}_{loadrundate}.rds")) # == Regressiontable
## Format estimate contrasts
Table[, 6:13] <- apply(Table[, 6:13], 2, pasteFormatted, t = Table$type)
## Format effects
Table[, 14:15] <- apply(Table[, 14:15], 2,
                        pasteFormatted, t = replace(as.character(Table$type), Table$type == "estimate", "percent"))

# Regression table --------------------------------------------------------
Regtable <- Table %>%
  filter(env != "P") %>%
  # mutate(opposite = ifelse(opposite, "!", "")) %>%
  mutate(n_clusters_ES = replace(n_clusters_ES, is.na(n_clusters_ES), ".")) %>%
  
  # keep only first row for factors 
  group_by(env, inventory, tax) %>%
  mutate_at(c("tax", "n_clusters_ES", "direction", "yr_betweenAvg", "yr_max"), keepFirstRow) %>% # not years!
  
  group_by(env, inventory) %>%
  mutate(environment = env) %>% # for later sorting
  mutate_at(c("env", "inventory"), keepFirstRow) %>%
  ungroup() %>%
  
  dplyr::select(environment, e = env, i = inventory, S = tax,
                adult_c = "(Intercept)_cond", "time:adult_c" = time_cond, juvenile_c = ontogenysmall_cond, "time:juvenile_c" = "time:ontogenysmall_cond",
                adult = "(Intercept)_disp", "time:adult" = time_disp, juvenile = ontogenysmall_disp, "time:juvenile" = "time:ontogenysmall_disp",
                "temp. shift" = timeeffect, "juv. divergence" = ontogenyeffect, direction, "no. presences" = n_pres,
                "years betw. avg." = yr_betweenAvg, "tot. period" = years, "no. clusters" = n_clusters_ES)


# dput(names(Regtable))
# digits <- c("e" = 0, "i" = 0, "Species" = 0,
#             "Intercept_c" = 3, "time_c" = 3, "juv._c" = 3, "time:juv._c" = 3,
#             "Intercept_d" = 3, "time" = 3, "juv." = 3, "time:juv." = 3,
#             "temp. shift" = 3, "ont shift" = 3, "o" = 0, "n_obs" = 0, "n_cl" = 0)

alignment <- c("l", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", "c", 
               "l", "r", "r", "r", "r")
# striped <- rep(seq(from = 0, to = nrow(Table)-5, by = 6), each = 3) + 1:3
opencurly <- "{"
closecurly <- "}"
caption <- "\\label{opencurly}tab:regtable_{str_replace(envname, fixed(' '), '_')}{closecurly} Regression estimates for {envname} with effects on mean $\\mu$ and dispersion parameter (precision) $\\phi$. The estimate for adults corresponds to the intercept, while other estimates are given as default treatment contrasts. For the models fitted to Spanish data, which had an additional fixed effect fitted for geographical clusters, contrast-equivalent marginal means and trends, averaged over the geographical clusters, are provided. In addition, marginal occurrence shifts on the response scale are given in percent of the sampled European gradient (shifts). The juvenile divergence is approximately equivalent to the following addition of model contrasts: \\code{opencurly}plogis(ad + juv) - plogis(ad){closecurly}. The temporal shift is approximately equivalent to \\code{opencurly}plogis(ad + juv +  time:ad + time:juv) - plogis(ad + juv){closecurly}. Stars * mark significant effects ($Pr(>|z|) = p < 0.05$). Arrows mark significant juvenile divergences that are either significantly opposite ($\\rightleftarrows$) or are significant in the same direction ($\\rightrightarrows$) as the temporal shift (one-sided test, $Pr(>|z|) = p < 0.05$). In addition the number of presences per life stage is given (juveniles:adults, number of total assessed plots in Supplementary Table \\ref{opencurly}tab:scope{closecurly}). Years between averages detail the difference between the average times of the presences at the second and the first sampling, while the total period spans all years with presences within country and species. Finally, the number of geographical clusters with presences in Spain is given."

Regtable_gdd <- dplyr::filter(Regtable, environment == "heat sum")
Regtable_cwb <- dplyr::filter(Regtable, environment == "water availability")

# Writes to files and returns html to viewer
catTable <- function(Tab, writetopaper = F){
  envname <- tolower(first(Tab$e))
  Tab <- dplyr::select(Tab, -e, -environment)
  
  n_rowspercountry <- c(table(filter(Table, tolower(env) == envname)$inventory))
  
  Tab_latex <- Tab
  names(Tab_latex) <- linebreak(names(Tab)) # converts \n linebreaks to latex witin table cell linebreaks
  
  
  Tab_latex %<>%
    dplyr::select(-i) %>%
    kable("latex", booktabs = T, longtable = T, caption = glue(caption), align = alignment, linesep = "") %>%
    add_header_above(c(" " = 1, "mean model estimates" = 4, "disp. model estimates" = 4, "shifts" = 3, " " = 4)) %>%
    kable_styling(font_size = 5, latex_options = c("repeat_header"), repeat_header_text = "\\textit{(continued)}", repeat_header_method = "replace") %>%
    # kable_styling(latex_options = "striped", stripe_index = striped, stripe_color = backgroundcolor) %>%
    pack_rows(index= c("Germany" = n_rowspercountry["Germany"], "Spain" = n_rowspercountry["Spain"], "Sweden" = n_rowspercountry["Sweden"])) %>%
    row_spec(0, angle = 80) %>%
    column_spec(1, italic = T)
  
  Tab_html <- kable(Tab, "html", booktabs = T, align = alignment, linesep = "") %>%
    add_header_above(c(" " = 2, "mean model estimates" = 4, "disp. model estimates" = 4, "shifts" = 3, " " = 4))
  
  
  cat(Tab_latex, file = glue("Publishing/Tables/{loadrunid}_Regression_Table_{str_replace(envname, fixed(' '), '_')}_{modelname}_{loadrundate}.tex"))
  cat(Tab_html, file = glue("Publishing/Tables/{loadrunid}_Regression_Table__{str_replace(envname, fixed(' '), '_')}_{modelname}_{loadrundate}.html"))
  if(writetopaper) cat(Tab_latex, file = glue("~/Documents/Studium/Projects/Papers/Mismatch/Tables/Regtable_{str_replace(envname, fixed(' '), '_')}.tex"))

  return(Tab_html)
}

catTable(Regtable_gdd, writetopaper = F)
# catTable(Regtable_gdd, writetopaper = T)

catTable(Regtable_cwb, writetopaper = F)
# catTable(Regtable_cwb, writetopaper = T)


# Scope table --------------------------------------------------------
Scopetable <- Table %>%
  droplevels() %>%
  dplyr::group_by(inventory) %>% # yr_betweenAvg
  summarize("selected species" = n_distinct(tax), "total period" = glue("{min(year1)}–{max(year2)}")) %>%
  dplyr::select("country" = "inventory", "selected species", "total period")

# Extract total number of assessed plots, not only presences of the chosen species
Plots_complete <- readRDS(glue("{inventoriespath}/Combined data/Env.rds")) %>%
  as.data.frame() %>% # drop geometry
  group_by(inventory) %>%
  summarize("assessed plots" = n_distinct(plotid)) %>%
  dplyr::select("assessed plots", "inventory")

i <- c(SE_NFI = "Sweden", DE_BWI = "Germany", ES_IFN = "Spain")

Scopetable <- cbind(Scopetable[match(i, Scopetable$country), ], Plots_complete[match(names(i), Plots_complete$inventory), "assessed plots"]) %>%
  dplyr::select("country", "assessed plots", "selected species", "total period")
  

Scopetab_latex <- Scopetable %>%
  kable("latex", booktabs = T, caption = glue("Sampling scope among the three NFI countries. Numbers of assessed plots and selected species show the intersecting subset that has been selected from different numbers of species and plots per inventory repetition. \\label{opencurly}tab:scope{closecurly}")) # %>%
  # add_header_above(c(" " = 3, "years between surveys" = 3))
  # kable_styling(font_size = 6)

Scopetab_html <- kable(Scopetable, "html") # %>%
  # add_header_above(c(" " = 3, "years between surveys" = 3))

cat(Scopetab_latex, file = glue("Publishing/Tables/{loadrunid}_Scope_Table_{modelname}_{loadrundate}.tex"))
# cat(Scopetab_latex, file = glue("~/Documents/Studium/Projects/Papers/Mismatch/Tables/Scopetable_new.tex"))

cat(Scopetab_html, file = glue("Publishing/Tables/{loadrunid}_Scope_Table_{modelname}_{loadrundate}.html"))


