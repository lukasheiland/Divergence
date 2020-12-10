
# Library  ------------------------------------------------------------------
library(extrafont)
# fonts()
# embed_fonts() ## Call file name after plotting to pdf to embed.
library(scales)
library(ggthemes)


# Colors ------------------------------------------------------------------

defaultcolor <- c("#13223F", "#FF4D26", "#5FEF8B", 
                  "#00C3FF", "#FCDA00", "#FF9050")


nicecolor <- setNames(defaultcolor, c("darkblue", "mennig", "mint", "lightblue", "yellow", "orange"))
ontogenycolor <- setNames(defaultcolor, c("big", "small", "big2", "small2", "big3", "small3"))
ontcolor <- setNames(defaultcolor, c("A", "J", "A2", "J2", "A3", "J3"))

options(ggplot2.discrete.color = defaultcolor)
options(ggplot2.discrete.fill = defaultcolor)

options(ggplot2.continuous.colour="viridis")
options(ggplot2.continuous.fill = "viridis")

## Set geom default color
# geomparam <- ls(pattern = '^geom_', envir = as.environment('package:ggplot2'))
# geomparam <- geomparam[!geomparam %in% c("geom_bin2d", "geom_count", "geom_freqpoly", "geom_histogram", "geom_jitter", "geom_qq", "geom_qq_line", 'geom_sf_label', 'geom_sf_text')]
basecolor <- defaultcolor[1]
geomparam <- c("geom_abline", "geom_area", "geom_bar", "geom_blank", "geom_boxplot", "geom_col", "geom_contour", "geom_contour_filled", "geom_crossbar", "geom_curve", "geom_density", "geom_density_2d", "geom_density_2d_filled", "geom_density2d", "geom_density2d_filled", "geom_dotplot", "geom_errorbar", "geom_errorbarh", "geom_function", "geom_hex", "geom_hline", "geom_label", "geom_line", "geom_linerange", "geom_map", "geom_path", "geom_point", "geom_pointrange", "geom_polygon", "geom_quantile", "geom_raster", "geom_rect", "geom_ribbon", "geom_rug", "geom_segment", "geom_sf", "geom_smooth", "geom_spoke", "geom_step", "geom_text", "geom_tile", "geom_violin", "geom_vline")
geomname <- gsub("geom_", "", geomparam)
lapply(geomname, update_geom_defaults, list(fill = basecolor, colour = basecolor))


# Themes ------------------------------------------------------------------

theme_ockham <- function(basetextsize = 12,
                   basefamily = ".SF Compact Text", # Georgia
                   baselwd = 1,
                   basecolor = "#13223F",
                   basealpha = 1,
                   axes = c("ranges", "axes", "none"),
                   suppressaxis = c("none", "x", "y"),
                   rangeframemap = NULL, # Provide deviating mapping for rangeframe.
                   legend = c("none", "left", "bottom", "right", "top"),
                   xlabangle = 0,
                   ...){
  
  defaultcolor <- c("#13223F", "#FF4D26", "#5FEF8B", 
                    "#00C3FF", "#FCDA00", "#FF9050")
  
  ## Set options
  options(ggplot2.discrete.color = defaultcolor)
  options(ggplot2.discrete.fill = defaultcolor)
  
  options(ggplot2.continuous.colour="viridis")
  options(ggplot2.continuous.fill = "viridis")
  
  ## Set geom defaults
  geomparam <- c("geom_abline", "geom_area", "geom_bar", "geom_blank", "geom_boxplot", "geom_col", "geom_contour", "geom_contour_filled", "geom_crossbar", "geom_curve", "geom_density", "geom_density_2d", "geom_density_2d_filled", "geom_density2d", "geom_density2d_filled", "geom_dotplot", "geom_errorbar", "geom_errorbarh", "geom_function", "geom_hex", "geom_hline", "geom_label", "geom_line", "geom_linerange", "geom_map", "geom_path", "geom_point", "geom_pointrange", "geom_polygon", "geom_quantile", "geom_raster", "geom_rect", "geom_ribbon", "geom_rug", "geom_segment", "geom_sf", "geom_smooth", "geom_spoke", "geom_step", "geom_text", "geom_tile", "geom_violin", "geom_vline")
  geomname <- gsub("geom_", "", geomparam)
  lapply(geomname, update_geom_defaults, list(fill = basecolor, colour = basecolor))
  
  if (!is.numeric(legend)) legend <- match.arg(legend)
  if (xlabangle > 5) xhjust <- 1 else xhjust <- NULL
  
  .theme <- theme(plot.margin = unit(c(t = 18, r = 16, b = 10, l = 10), "pt"),
                  plot.background = element_blank(),
                  strip.background = element_blank(),
                  
                  panel.background = element_blank(),
                  panel.border = element_blank(),
                  panel.grid = element_blank(),
                  
                  text = element_text(size = basetextsize,
                                      family = basefamily),
                  title = element_text(family = basefamily, colour = basecolor), # all titles
                  plot.title = element_text(hjust = 0.5), # main title
                  strip.text = element_text(family = basefamily, colour = basecolor), # facet title
                  strip.text.x = element_text(margin = margin(t = 3, r = 1, b = 4, l = 1)),

                  axis.line = element_blank(), # Use geom_rangeframe instead.
                  axis.text = element_blank(),
                  axis.title = element_blank(),
                  axis.ticks = element_blank(),
                  
                  legend.background = element_blank(),
                  legend.box.background = element_rect(fill = NULL, color = scales::alpha(basecolor, basealpha)),
                  legend.position = legend,
                  legend.key = element_blank(),
                  ...
          )
  
  if (xlabangle != 0){
    .theme <- .theme + theme(axis.text.x = element_text(angle = xlabangle, hjust = xhjust))
    }
  
  if (axes[1] != "none") {
    .theme <- .theme +
      theme(axis.ticks = element_line(color = scales::alpha(basecolor, basealpha), size = baselwd*0.55),
            axis.text = element_text(color = scales::alpha(basecolor, basealpha), size = basetextsize*0.8),
            axis.text.x = element_text(margin = margin(t = 2.2, b = 1, unit = "pt")),
            axis.text.y = element_text(margin = margin(r = 2, l = 1, unit = "pt")),
            axis.title = element_text(color = basecolor, size = basetextsize),
            axis.title.x = element_text(margin = margin(t = 7, b = 1, unit = "pt")),
            axis.title.y = element_text(margin = margin(r = 7, l = 1, unit = "pt"))
            )
    
    if(axes[1] == "ranges"){
      .theme <- list(.theme,
                     geom_rangeframe(mapping = rangeframemap, color = basecolor, alpha = basealpha, size = baselwd)
                     )
      }
    
    if(axes[1] == "axes"){
      .theme <- list(.theme + theme(axis.line = element_line(color = scales::alpha(basecolor, basealpha), size = baselwd)) #,
                     # scale_y_continuous(expand = c(0, 0))
      )
      }
    
    # if(suppressaxis[1] == "y") {
    #   .theme <- .theme +
    #     theme(axis.line.y = element_blank(), # Use geom_rangeframe instead.
    #           axis.text.y = element_blank(),
    #           # axis.title.y = element_blank(),
    #           axis.ticks.y = element_blank()
    #     )
    #   }
    # 
    # if(suppressaxis[1] == "x") {
    #   .theme <- .theme +
    #     theme(axis.line.x = element_blank(), # Use geom_rangeframe instead.
    #           axis.text.x = element_blank(),
    #           # axis.title.y = element_blank(),
    #           axis.ticks.x = element_blank()
    #     )
    #   }
  }
  
  return(.theme )
}


theme_empty <- function(){
  theme(panel.grid = element_line(colour = "transparent"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
  }

