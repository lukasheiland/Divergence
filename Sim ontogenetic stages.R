# The Library ---------------------------------------------------------
library(deSolve)
library(extraDistr)
library(caTools) # fastest and feature complete running mean
library(scales)

library(tidyr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(glue)

library(glmmTMB)
library(emmeans)
library(ggbeeswarm)
library(rescale)

library(kableExtra)


# Set working directory ---------------------------------------------------
setwd(rprojroot::find_root(rprojroot::criteria$is_git_root))


# Source ------------------------------------------------------------------
if (file.exists("Theme.R")) source("Theme.R") else theme_ockham <- function(rangeframemap) NULL


# Set up the environmental gradient ------------------------------------
n_env <- 300
environment <- seq(0, 1, length.out = n_env)


# Stage simulation model -------------------------------------------------------------

## "Gaussian" environmental filter with y = 1 at optimum.
filtering <- function(env, optimum, width) exp(-0.5 * ((env - optimum) / width)^2)

## "Quadratic" environmental upside down filter.
# filteringM <- function(env, optimum, width) 1/width * ((env - optimum))^2

## "Gaussian" environmental upside down filter with 0 at minimum.
filteringM <- function(env, optimum, width) 1 - exp(-0.5 * ((env - optimum) / width)^2)

#### ODE systems for two logistically behaving stages
## r: regeneration rate of juveniles.
## trans: transitioning fraction from J to A [formerly b]
## s: survival rates of adults.
## K_*: Carrying capacity at the environment optimum. Mortality caused by environment is implicitly included in this.
## k_*: vector of carrying capacitites k[1:n_env]
## J/A: states of pops
## D: Diaspore stage state. Basically A, spread out with ...
## d: rolling average by parameter for D.



### 1. Logistic increase and reduction model
## (i.e. logistic death rates in adults and birth rates in juveniles.
## Transition is independent of K, and thus independent of env.

odeLogistic <- function(t, state, par, n_env, filterfun_j, filterfun_a) {
  J <- state[1:n_env]
  A <- state[(n_env+1):(2*n_env)]
  
  with(par,
       {
         k_j <- K_j * filterfun_j(t) # Vector of carrying capacities, k(env).
         k_a <- K_a * filterfun_a(t) #
         
         ## This fun idea for dispersal wouldn't work to be integrated.
         ## Gaussian dispersal kernel. Where the number of Diaspores D are equal to the number of adults A in an adjacent cell of env.
         # i_d <- round(rtnorm(n_env, 1:n_env, sd = disp, a = 1, b = n_env))
         # D <- A[i_d] # a < x <= b!s
         
         ## Dispersal is rolling window average of A, center-aligned. + constant rate seedbank
         D <- runmean(A, d, endrule = "mean", align = "center") # + seedbank
         
         dJ <- r*D * (1 - J / k_j) - trans*J ## Temporal change of juveniles J.
         dA <- trans*J - m*A * (A / k_a) ## Temporal change of adults A.

         return(list(c(dJ, dA)))
       }
  )
}



### 2. Two stage model with a 1/Normal mortality "filter".
## Here, K is constant over environment, limiting growth rates just as an expression of density dependence.
## But mortality, m, is a function of env.

# ## Here the scaling of the filtering is through division, leading to an inverted bathtub with infinite mortality at the margins.
# odeMortality <- function(t, state, par, n_env, filterfun_j, filterfun_a) {
#   J <- state[1:n_env]
#   A <- state[(n_env+1):(2*n_env)]
#   
#   with(par,
#        {
#          m_j <- M_j/filterfun_j(t) # filterfun is 1 at the maximum.
#          m_a <- M_a /filterfun_a(t) #
#     
#          D <- runmean(A, d, endrule = "mean", align = "center") + seedbank
#          dJ <- (r*D - m_j*J) * (1 - J/K_j) - trans*J ## Temporal change of juveniles J.
#          dA <- trans*J - m_a*A * (A/K_a) ## Temporal change of adults A.
#          
#          return(list(c(dJ, dA)))
#        }
#   )
# }




### 3. Two stage model with a normal mortality "filter".
## Here, K is constant over environment, limiting growth rates just as an expression of density dependence.
## But mortality, m, is a function of env.

odeMortality <- function(t, state, par, n_env, filterfun_j, filterfun_a) {
  J <- state[1:n_env]
  A <- state[(n_env+1):(2*n_env)]
  
  with(par,
       {
         m_j <- M_j * filterfun_j(t) # filterfun is 1 at the maximum. Was like this before for Gaussian: M_j * (1 - filterfun_j(t))
         m_a <- M_a * filterfun_a(t) #
         
         k_j <- K_j * (1 - filterfun_j(t))
         k_a <- K_a * (1 - filterfun_a(t))
         
         D <- runmean(A, d, endrule = "mean", align = "center")
         dJ <- (r*D - m_j*J- trans*J) * (1 - J/k_j) ## Temporal change of juveniles J.
         dA <- (trans*J - m_a*A) * (1 - A/k_a) ## Temporal change of adults A.
         
         return(list(c(dJ, dA)))
       }
  )
}

## Same model, but the stages are independent of each other.
odeMortalityNiche <- function(t, state, par, n_env, filterfun_j, filterfun_a) {
  J <- state[1:n_env]
  A <- state[(n_env+1):(2*n_env)]
  
  with(par,
       {
         m_j <- M_j * filterfun_j(t) # filterfun is 1 at the maximum. Was like this before for Gaussian: M_j * (1 - filterfun_j(t))
         m_a <- M_a * filterfun_a(t) #
         
         k_j <- K_j * (1 - filterfun_j(t))
         k_a <- K_a * (1 - filterfun_a(t))
         
         dJ <- (r*J - m_j*J) * (1 - J/k_j) ## Temporal change of juveniles J.
         dA <- (r*A - m_a*A) * (1 - A/k_a) ## Temporal change of adults A. ##!!! used r instead of trans because otherwise reproduction would be to small.
         
         return(list(c(dJ, dA)))
       }
  )
}



### X. Competition model
# odeLogistic_comp <- function(t, y, parms)
# {
#   with(append(y, parms),
#        {
#          ### new parameters:
#          ## c_* competition of
#          ## d_* death rate [formerly d, only for A: d_a]
#          
#          dJ <- r*A*(1 - (c_a*A + c_j*J) / k_j) - trans*J
#            ## Temporal change of juveniles J.
#          dA <- trans*J - (d_a*A*(c_a*A + c_j*J)) / k_a
#            ## Temporal change of adults A.
#          return(list(c(dJ, dA)))
#        }
#   )
# }




# Integration of the ODE model --------------------------------------------
## ... proviing changing filters per time step.

simulateAbundance <- function(model,
                              pars,
                              times,
                              env = environment,
                              state = c(
                                J = rep_len(0.1, length(env)),
                                A = rep_len(0.1, length(env))
                              ),
                              
                              ## The following parameters inform the filter functions filterfun_*
                              filter_j = c(m = 0.5, s = 0.13),
                              filter_a = c(m = 0.5, s = 0.13),
                              filterslope = c(m = 0, s = 0),
                              
                              returnstate = F ## Whether to return a state vector of the last time in times, instead of data.frame
                              ) {
  
  filtering <- if(identical(model, odeMortality)) filteringM else filtering
    
  filterfun_j <- function(t) filtering(env, filter_j[1] + filterslope[1] * t, filter_j[2] + filterslope[2] * t)
  filterfun_a <- function(t) filtering(env, filter_a[1] + filterslope[1] * t, filter_a[2] + filterslope[2] * t)
  
  Out <- ode(y = state,
             times = c(0, times),
             func = model,
             parms = pars,
             n_env = length(env),
             filterfun_j = filterfun_j,
             filterfun_a = filterfun_a)
  
  if(returnstate) return(Out[length(times) + 1, -1])
  

  A <- tidyr::pivot_longer(as.data.frame(Out),
                           starts_with("J") | starts_with("A"),
                           names_to = "ontogeny",
                           values_to = "abundance") %>%
    mutate(ontogeny = substr(ontogeny, 1, 1)) %>%
    filter(time != 0) %>% 
    bind_cols(env = rep(env, times = 2*length(times))) %>%
    mutate(filter = ifelse(ontogeny == "J", ## just the values of the filter functions
                           filterfun_j(time),
                           filterfun_a(time)),
           K = ifelse(ontogeny == "J", ## K dependent on filter 
                  filter * pars$K_j,
                  filter * pars$K_a),
           M = ifelse(ontogeny == "J", ## M dependent on filter
                      (1 - filter) * max(abundance), # Truly: (1 - filter) * pars$M_a
                      (1 - filter) * max(abundance))
           ) %>%
    group_by(time, env) %>%
    mutate(total = sum(abundance)) %>%
    ungroup() %>%
    group_by(time) %>%
    mutate(lower_total = quantile(sample(env, size = 1000*length(env), replace = T, prob = replace(total, total < 0, 0)), 0.025), upper_total = quantile(sample(env, size = 1000*length(env), replace = T, prob = replace(total, total < 0, 0)), 0.975)) %>%
    ungroup() %>%
    mutate(total_last = rep(total[time == max(time)], length(times)))
  
  if(identical(model, odeMortality)) {
    nichemodel <- odeMortalityNiche
    Niche <- ode(y = c(J = rep_len(0.1, length(env)),
                       A = rep_len(0.1, length(env))
                       ),
                 times = c(0, times), # same times as above
                 func = nichemodel,
                 parms = pars,
                 n_env = length(env),
                 filterfun_j = filterfun_j,
                 filterfun_a = filterfun_a)
    
    A_niche <- tidyr::pivot_longer(as.data.frame(Niche),
                                   starts_with("J") | starts_with("A"),
                                   names_to = "ontogeny",
                                   values_to = "abundance") %>%
      mutate(ontogeny = substr(ontogeny, 1, 1)) %>%
      dplyr::filter(time != 0) %>%
      dplyr::select(niche = abundance)
    
    A <- cbind(A, A_niche)
  }
  
  attr(A, "simparameters") <- pars
  attr(A, "filterparameters") <- c(m_j = filter_j[[1]], s_j = filter_j[[2]], m_a = filter_a[[1]], s_a = filter_a[[2]], b_m = filterslope[[1]], b_s = filterslope[[2]])
  return(A)
}


# Data summary -----------------------------------------------------

## Transforms abundance data into presence absence data through resampling.
getPresData <- function(A, env = environment, sampletimes = 8) {
  A %<>%
    filter(abundance > 0) %>%
    filter(env > 0 & env < 1) %>%
    slice_sample( n = sampletimes*length(env), weight_by = abundance, replace = T) %>%
    mutate(time = scales::rescale(time, to = c(-0.5, 0.5)))
  return(A)
}


## Fit beta model to presence data
getFit <- function(P) glmmTMB(env ~ time * ontogeny, dispformula = ~ time * ontogeny, family = beta_family, data = P)



## Calculate summary statistics
getStats <- function(A, env = environment) {
  J <- filter(A, ontogeny == "J" & time == max(time))
  A <- filter(A, ontogeny == "A" & time == max(time))
  
  prob_j <- J$abundance / sum(J$abundance)
  prob_a <- A$abundance / sum(A$abundance)
  prob_j[prob_j < 0] <- 0
  prob_a[prob_a < 0] <- 0
  
  # sample
  env_j <- sample(env, size = 1000*length(env), replace = T, prob = prob_j)
  env_a <- sample(env, size = 1000*length(env), replace = T, prob = prob_a)
  
  range_j <- max(env[J$abundance > 0.1]) - min(env[J$abundance > 0.1])
  range_a <- max(env[A$abundance > 0.1]) - min(env[A$abundance > 0.1])
  
  upper_j <- as.numeric(quantile(env_j, 0.95))
  median_j <- as.numeric(quantile(env_j, 0.5))
  lower_j <- as.numeric(quantile(env_j, 0.05))
  upper_a <- as.numeric(quantile(env_a, 0.95))
  median_a <- as.numeric(quantile(env_a, 0.5))
  lower_a <- as.numeric(quantile(env_a, 0.05))
  
  Stat <- data.frame(ontogeny = c("J", "A"),
                     mean = c(mean(env_j), mean(env_a)),
                     sd = c(sd(env_j), sd(env_a)),
                     range = c(range_j, range_a),
                     inr = c(upper_j - lower_j, upper_a - lower_a),
                     upper = c(upper_j, upper_a),
                     median = c(median_j, median_a),
                     lower = c(lower_j, lower_a),
                     
                     stringsAsFactors = F
                     )
  return(Stat)
}


## Returns a data.frame of the ontogeny effect:
## Same as getOntogenyEffect() in analyses but with less overhead.
## (small intercept - big intercept) [hence: revpairwise!] by inventory,
## The null whether == 0 is tested (disp model on log scale, which means on response scale == 1) .
getSimOntogenyEffect <- function(fit,
                                 glmmTMBcomponent = "cond", # Which model component to use, dispersion or "conditional" (µ)? Passed on to https://github.com/glmmTMB/glmmTMB/blob/78856c092e286f4d449af707bd7fdfaf638957c3/glmmTMB/R/emmeans.R Alternatives: "cond", "disp".
                                 attime = NULL, # c(-1, 0, 1),
                                 abs = FALSE,
                                 transform = "response", # alternative: "none"
                                 averageoverinventories = FALSE) {
  
  glmmTMBcomponent <- glmmTMBcomponent ## https://stackoverflow.com/questions/18136720/why-missing-and-default-arguments-are-not-working-in-functions-called-by-lapp
  
  if (is.null(attime)) {
    ont <- emmeans::emmeans(fit,
                            specs = revpairwise ~ ontogeny,
                            transform = transform,
                            component = glmmTMBcomponent)$contrasts
    
  } else {
    ont <- emmeans::emmeans(fit,
                            specs = revpairwise ~ ontogeny | inventory:time,
                            at = list(time = attime),
                            transform = transform,
                            component = glmmTMBcomponent
    )$contrasts
  }
  
  CL <- as.data.frame(confint(ont))[c("lower.CL", "upper.CL")] # level 0.95;
  Out <- cbind(as.data.frame(ont), CL) %>% 
    mutate(signif = p.value < 0.05) %>%
    mutate(star = symnum(p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.05, 1), symbols = c("*", "n.s."))) %>%
    cbind(effect = "ontogeny")
  
  if (abs & !glmmTMBcomponent == "disp") {
    inverted <- Out$estimate < 0
    lower <- Out$lower.CL
    Out$lower.CL[inverted] <- Out$upper.CL[inverted]
    Out$upper.CL[inverted] <- lower[inverted]
    Out$estimate <- abs(Out$estimate)
  }
  
  return(Out)
}


## Returns a data.frame including the slopes (time effects) by ontogeny and inventory,
## Same as getTimeEffect() in analyses but with less overhead.
## The null whether == 0 is tested (disp model on log scale, which means on target scale == 1) .
getSimTimeEffect <- function(fit,
                             glmmTMBcomponent = "cond", # Which model component to use, dispersion or "conditional" (µ)? Passed on to https://github.com/glmmTMB/glmmTMB/blob/78856c092e286f4d449af707bd7fdfaf638957c3/glmmTMB/R/emmeans.R Alternatives: "cond", "disp".
                             ontogenylevel = c("J", "A"),
                             abs = FALSE,
                             transform = "response", # alternative: "none"
                             averageoverinventories = FALSE) {
  
  glmmTMBcomponent <- glmmTMBcomponent ## https://stackoverflow.com/questions/18136720/why-missing-and-default-arguments-are-not-working-in-functions-called-by-lapp
  
  timetrend <- emmeans::emtrends(fit,
                                 specs = ~ ontogeny,
                                 var = "time",
                                 transform = transform,
                                 component = glmmTMBcomponent)
  
  timetest <- emmeans::test(timetrend,
                            null = 0,
                            joint = F,
                            by = c("ontogeny")) # "by" ensures independent tests!
  
  
  
  CL <- as.data.frame(timetrend)[c("lower.CL", "upper.CL")] # level 0.95
  Out <- cbind(as.data.frame(timetest), CL)
  
  if (abs & !glmmTMBcomponent == "disp") {
    inverted <- Out$time.trend < 0
    lower <- Out$lower.CL
    Out$lower.CL[inverted] <- Out$upper.CL[inverted]
    Out$upper.CL[inverted] <- lower[inverted]
    Out$time.trend <- abs(Out$time.trend)
  }
  
  Out <- dplyr::filter(Out, ontogeny %in% ontogenylevel) %>%
    dplyr::rename(estimate = "time.trend") %>%
    mutate(signif = p.value < 0.05) %>%
    mutate(star = symnum(p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.05, 1), symbols = c("*", "n.s."))) %>%
    cbind(effect = "time")
  
  return(Out)
}





# Plotting ----------------------------------------------------------------

## Plot lines from well defined data.frames
plotLines <- function(Abundance, Presence, Effects, fit, spread = 0.3) {
  
  A <- Abundance %>%
    mutate(time = scales::rescale(time, to = c(-0.5, 0.5))) %>%
    mutate(abundance = if_else(
      time == min(time),
      scales::rescale(abundance, to = c(min(time), min(time) + spread)),
      scales::rescale(abundance, to = c(max(time), max(time) + spread))
    ))
  
  lshiftstart <- plogis(fit$fit$par[1])
  lshiftend <- plogis(fit$fit$par[1] + fit$fit$par[3])
  
  tshiftstart <- predict(fit, newdata = list(time = -0.5, ontogeny = "J"), type = "response")
  tshiftend <- predict(fit, newdata = list(time = 0.5, ontogeny = "J"), type = "response")
  
  lshift <- Effects[Effects$effect == "ontogeny", "estimate"]
  tshift <- Effects[Effects$effect == "time" & Effects$ontogeny == "J", "estimate"]
  
  lstar <- Effects[Effects$effect == "ontogeny", "star"]
  tstar <- Effects[Effects$effect == "time" & Effects$ontogeny == "J", "star"]
  
  llabel <- glue("life stage shift
                 {sprintf('%.1f%s', lshift * 100, '%')} {lstar}")
  tlabel <-  glue("temporal shift
                  {sprintf('%.1f%s', tshift * 100, '%')} {tstar}")
  
  
  
  ggplot(Presence, aes(y = env, x = time, col = ontogeny, group = ontogeny, fill = ontogeny)) +
    theme_ockham() +
    ylim(min(Presence$env) + 0.1, max(Presence$env) - 0.1) +
    
    ## arrow for life stage shift
    geom_segment(aes(x = 0, xend = 0, y = lshiftstart, yend = lshiftend), arrow = arrow(angle = 20, ends = "last", length = unit(6, "pt")), size = 0.1, col = 1) +
    geom_curve(aes(x = 0, y = lshiftstart + 0.03, xend = 0, yend = lshiftstart + 0.15), curvature = -0.4, angle = 60, lwd = 0.004, col = 1) +
    annotate("label", x = 0, y = lshiftstart + 0.15, label = llabel, fill = 'white', label.size = NA) +
    
    geom_segment(aes(x = -0.48, xend = 0.48, y = tshiftstart, yend = tshiftstart), lty = 3, size = 0.08, col = 1) +
    geom_segment(aes(x = 0.49, xend = 0.49, y = tshiftstart, yend = tshiftend), arrow = arrow(angle = 20, ends = "last", length = unit(6, "pt")), size = 0.1, col = 1) +
    geom_curve(aes(x = 0.49, y = tshiftstart - 0.01, xend = 0.2, yend = tshiftstart - 0.15), curvature = 0.4, angle = 120, lwd = 0.004, col = 1) +
    annotate("label", x = 0.1, y = tshiftstart - 0.15, label = tlabel, fill = 'white', label.size = NA) +
    
    
    geom_smooth(method = "lm", fill = NA) + # lwd = 2
    geom_quasirandom(alpha = 0.15, width = spread, groupOnX = T, dodge.width = 0.03, size = 0.01) + # with dodge.width != 0 groups will be separated
    geom_path(data = A, aes(y = env, x = abundance, group = interaction(time, ontogeny), color = ontogeny))
  
  # geom_point(position = position_jitter(width = 0.05), alpha = 0.01)
}


## Plot a well-defined Abundance data.frame
plotAbundance <- function(A, title = NULL, subtitle = NULL, plotfilter = "K", plotrange = T, plotcenter = T, plottotal = F, plotyaxis = F, themefun = theme_ockham, axes = "ranges", ...){
  
  A <- left_join(A, getStats(A), by = "ontogeny")
  fewtimes <- length(unique(A$time)) < 4
  if(fewtimes) A %<>% mutate(arrowypos = -1 - (as.numeric(as.factor(A$time))+1))
  A %<>% mutate(medianxpos = max(abundance))

  plot <- ggplot(A, aes(env, abundance, color = ontogeny, group = ontogeny)) +
    scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3, digits = 0)) + # labels = scales::scientific
    xlab("environment at t = 1") +
    { if (any(is.character(title), is.character(subtitle))) ggtitle(title, subtitle) } +
    { if (length(unique(A$time)) > 1) aes(alpha = time, group = interaction(ontogeny, time)) } +
    {if (fewtimes) geom_path(size = 1.5) else geom_path()} +
    scale_alpha(range = c(0.3, 1)) +
    { if (is.character(plotfilter)) geom_path(aes_(~env, as.name(plotfilter)), linetype = 1, size = 0.7, color = "white") } + # A white background for the filter.
    { if (is.character(plotfilter)) geom_path(aes_(~env, as.name(plotfilter)), linetype = 3) } +
    { if (is.character(plotfilter) &  plotyaxis) themefun(rangeframemap = aes_(~env, as.name(plotfilter)), axes = axes, ...) else themefun(axes = axes, rangeframemap = aes(y = NULL), suppressaxis = "y", ...) } +
    
    ## Reduction functions like min/max/median are just hacks to reduce plotting resources.
    { if (plotcenter & fewtimes) geom_segment(aes(x = median, xend = median, y = medianxpos, yend = medianxpos*1.15, col = ontogeny, group = time), size = 0.4) } +
    { if (plotcenter & fewtimes) geom_label(aes(label = "50%", x = median, y = medianxpos*1.15, col = ontogeny, group = time), fill = "white", label.size = NA) } +
    
    { if (plotrange & !plottotal & fewtimes) geom_segment(aes(x = lower, xend = upper, y = -1, yend = -1), arrow = arrow(angle = 90, ends = "both", length = unit(6, "pt")), size = 0.4) } +
    { if (plotrange & !plottotal & fewtimes) geom_text(aes(label = "5%", x = lower, y = 0, col = ontogeny)) } +
    { if (plotrange & !plottotal & fewtimes) geom_text(aes(label = "95%", x = upper, y = 0, col = ontogeny)) } +
    
    { if (plotrange & plottotal & fewtimes) geom_segment(aes(x = lower_total, xend = upper_total, y = arrowypos, yend = arrowypos, group = time), arrow = arrow(angle = 90, ends = "both", length = unit(6, "pt")), size = 0.4, colour = "grey") } +
    { if (plotrange & plottotal & fewtimes) geom_text(aes(label = "5%", x = lower_total, y = arrowypos+1, group = time), colour = "grey") } +
    { if (plotrange & plottotal & fewtimes) geom_text(aes(label = "95%", x = upper_total, y = arrowypos+1, group = time), colour = "grey") } +
    
    { if (plottotal) geom_path(aes(env, total), linetype = "longdash", color = "grey")} +
    { if (plottotal & plotyaxis) themefun(rangeframemap = aes(env, total), axes = axes, ...) else themefun(axes = axes, rangeframemap = aes(y = NULL), suppressaxis = "y", ...) } +
    
    scale_color_manual(values = c(A = "#13223F", J = "#FF4D26",
                                  A2 = "#002166", J2 = "#EA2B00",
                                  A3 = "#3D4A66", J3 = "#FF8266")) +
    theme(aspect.ratio = 1/2)
  
  ## Attach "attributes" to the ggplot object, for later extraction
  plot$simparameters <- attr(A, "simparameters")
  plot$filterparameters <- attr(A, "filterparameters")
  plot$Abundance <- A

  return(plot)
}


#### Here are just wrapper functions for quick combined simulation and plotting:
plotOntogenetic <- function(model, pars, filter_j = c(0.5, 0.13), filter_a = c(0.5, 0.13), times = 10000,
                            title = "X)", plottotal = F, plotrange = F, axes = "ranges", plotregression = F, ...) {
  
  A <- simulateAbundance(model, pars, filter_j = filter_j, filter_a = filter_a, times = times)
  plotfilter <- if(identical(model, odeMortality)) "niche" else "K" # alternative "niche"
  plot <- plotAbundance(A, plottotal = plottotal, plotrange = plotrange, axes = axes, subtitle = title, plotfilter = plotfilter, ...)
  
  if(plotregression) {
    P <- getPresData(A)
    fit <- getFit(P)
    E <- bind_rows(getSimOntogenyEffect(fit),
                   getSimTimeEffect(fit))
    lineplot <- plotLines(A, P, E, fit)
    
    plot <- ggdraw(plot) +
      draw_plot(lineplot, x = 0.55, y = 0.36, width = 0.45, height = 0.65) # + draw_plot_label(c("A", "B"), c(0, 0.45), c(1, 0.95), size = 12)
  }
  
  ## Attach "attributes" to the ggplot object, for later extraction
  plot$simparameters <- pars
  plot$Abundance <- A
  plot$filterparameters <- attr(A, "filterparameters")
  
  return(plot)
}


plotTemporal <- function(model, pars, filter_j = c(0.5, 0.13), filter_a = c(0.5, 0.13), filterslope = c(0, 0), times,
                         title = "X)", plottotal = F, plotrange = T, axes = "ranges", plotregression = F, ...) {
  
  steadystate <- simulateAbundance(model, pars, filter_j = filter_j, filter_a = filter_a, times =  c(10000), returnstate = T)
  A <- simulateAbundance(model, pars, filter_j = filter_j, filter_a = filter_a, filterslope = filterslope, times = times, state = steadystate)
  plotfilter <- if(identical(model, odeMortality)) "niche" else "K" # alternative "niche"
  plot <- plotAbundance(A, plottotal = plottotal, plotrange = plotrange, axes = axes, subtitle = title, plotfilter = plotfilter, ...)
  
  if(plotregression) {
    P <- getPresData(A)
    fit <- getFit(P)
    E <- bind_rows(getSimOntogenyEffect(fit),
                   getSimTimeEffect(fit))
    lineplot <- plotLines(A, P, E, fit)
    
    plot <- ggdraw(plot) +
      draw_plot(lineplot, x = 0.55, y = 0.36, width = 0.45, height = 0.65) # + draw_plot_label(c("A", "B"), c(0, 0.45), c(1, 0.95), size = 12)
  }
  
  ## Attach "attributes" to the ggplot object, for later extraction
  plot$simparameters <- pars
  plot$Abundance <- A
  plot$filterparameters <- attr(A, "filterparameters")
  
  return(plot)
}
  


# Set up for simulation  ----------------------------------------------------------

#### Reminder of the defaults:
formals(simulateAbundance)

#### Set parameters.
p_log <- list(r = 0.2, # Regeneration, acting on J.
              trans = 0.008, # Fraction transitioning from J to A. # 0.0066
              m = 0.03,
              K_j = 20, # Total carrying capacity, will be stretched out over env.
              K_a = 10,
              d = 50 # d: integer, width of rolling window to apply to adults for dispersal.
              # seedbank = 0
              )

# p_mort_old <- within(p_log, {M_a <- 0.02 # minimal mortality at the optimum
#                              M_j <- 0.02,
#                              m = 0.03
# })

p_mort <- list(r = 0.2, # Regeneration, acting on J.
               trans = 0.02, # Fraction transitioning from J to A. # 0.0066
               K_j = 10000, # Total carrying capacity, will be stretched out over env.
               K_a = 7000,
               d = 50, # d: integer, n_cells width of rolling window to apply to adults for dispersal.
               M_j = 0.1, # Maximal mortality at the margins. The mortality niche ranges from 0 at the optimum to M. 
               M_a = 0.1) # Has to be bigger than r or trans, then it will steer pop towards zero.


t_long <- (1:30)^2
t_linear <- seq(10, 500, by = 10)
t_linear2 <- seq(600, 5000, by = 100)
t <- c(t_linear, t_linear2)



#############################################################################
# Simulate cases -----------------------------------------------------------#
#############################################################################

# ## Minimal example of complete tool chain, which is wrapped in plotOntogenetic(), and plotTemporal()
# t_fit <- c(100, 500)
# state_A_M_ont <- simulateAbundance(odeMortality,
#                                    within(p_mort, {trans <- 0.005; r <- 3;}), times = c(5000),
#                                    filter_j = c(0.45, 0.13), filter_a = c(0.5, 0.13),
#                                    returnstate = T)
# A_M_ont <- simulateAbundance(odeMortality, state = state_A_M_ont,
#                              within(p_mort, {trans <- 0.005; r <- 3; M_a <- 0.05}),
#                              filter_j = c(0.45, 0.13), filter_a = c(0.5, 0.13),
#                              filterslope = c(0.0001, 0), times = t_fit)
# mainplot <- plotAbundance(A_M_ont, plotfilter = "niche"); mainplot
# P_M_ont <- getPresData(A_M_ont)
# fit_M_ont <- glmmTMB(env ~ time * ontogeny, dispformula = ~ time * ontogeny, family = beta_family, data = P_M_ont) # see also: getFit()
# 
# E_M_ont <- bind_rows(getSimOntogenyEffect(fit_M_ont),
#                      getSimTimeEffect(fit_M_ont))
# regressionplot <- plotLines(A_M_ont, P_M_ont, E_M_ont, fit_M_ont); regressionplot
# 
# ggdraw(mainplot + theme_half_open(12)) +
#   draw_plot(regressionplot, .45, .45, .5, .5)



# I) Ontogenetic change -------------------------------------------------------

##### CASE 0) No change
## 
A_L_same <- simulateAbundance(odeLogistic, p_log,
                                times = t_linear)
plot_L_same <- plotAbundance(A_L_same, subtitle = "X) Equal carrying capacity", plotrange = F, plotcenter = F)
getStats(A_L_same)

##
plot_M_same <- plotOntogenetic(odeMortality, within(p_mort, {K_j <- K_a <- 10000}), filter_j = c(0.5, 0.12), filter_a = c(0.5, 0.12),
                               title =  "Ontogenetically equal mortality rates and response", plotcenter = F, plotrange = T)
# within(p_mort, {K_j <- 20; K_a <- 20})
## Assuming constant growth rates, J and A have the same mortality "niche". Still the resulting occurrence differs in width. Mortality niche depending on environment while being constant over life-history, will lead to relatively lower abundance at high mortality conditions in later life stages.


## Expanding occurrence ontogenetically with same K niche does not work.
# A_L_same_structJnarrower <- simulateAbundance(odeLogistic,
#                                               within(p_log, {trans <- 0.002; m = 0.0003; K_j <- 20; K_a <- 10}),
#                                            filter_j = c(0.5, 0.1), filter_a = c(0.5, 0.1),
#                                            times = t_linear)
# plot_L_same_structJnarrower <- plotAbundance(A_L_same_structJnarrower, subtitle = "X) Ontogenetically equal K niche", axes = T, plotrange = T)




##### CASE Ont1) Change of µ over ontogeny
##
plot_L_ont_m <- plotOntogenetic(odeLogistic, p_log, filter_j = c(0.4, 0.1), filter_a = c(0.6, 0.1),
                                title = "X) Ontogenetically shifting carrying capacity niche")
##
set.seed(100)
plot_M_ont_m <- plotOntogenetic(odeMortality, p_mort, filter_j = c(0.5, 0.15), filter_a = c(0.3, 0.15), times = c(1e4, 1e4 + 100),
                                title = "Ontogenetically shifting niche", plotcenter = F, plotregression = T)


##### CASE Ont2) Change of filter σ over ontogeny
##
plot_L_ont_s <- plotOntogenetic(odeLogistic, p_log, filter_j = c(0.5, 0.1), filter_a = c(0.5, 0.4),
                                title = "X) Ontogenetically expanding carrying capacity niche", plotcenter = F)

# ## J > A with K niche is possible
# plot_L_ont_s <- plotOntogenetic(odeLogistic, p_log, filter_j = c(0.5, 0.2), filter_a = c(0.5, 0.05),
#                                 title = "X) Ontogenetically contracting carrying capacity niche", plotrange = T)

## J < A with m niche is almost(?) impossible to produce (maybe with very small effect and extreme assumptions)
plot_M_ont_s <- plotOntogenetic(odeMortality, within(p_mort, {d <- 0}), filter_j = c(0.5, 0.12), filter_a = c(0.5, 0.15),
                                title = "X) Ontogenetically expanding mortality niche", plotcenter = F)
## Juvenile occurrence limiting adult occurrence, despite wider adult potential. Vice versa, equilibrium of the two niches.





##### CASE Ont3) Change of filter sigma and µ over ontogeny
# ## boooooring
# plot_L_ont_ms <- plotOntogenetic(odeLogistic, p_log,  filter_j = c(-1, 0.5), filter_a = c(0, 1), , plotrange = T)

plot_M_ont_ms <- plotOntogenetic(odeMortality, p_mort, filter_j = c(0.45, 0.12), filter_a = c(0.53, 0.15),
                                title = "X) Ontogenetically expanding and shifting mortality niche", plotrange = T)





# II) Temporal change -------------------------------------------------------

## Linear dependency of the "niche" function on time.
## filterslope[1] = effect on mean
## filterslope[2] = effect on spread
##

## Set of times
t_temp <- c(1,100)
t_fit <- c(100, 500) # This is

#### CASE Temp1) Temporal change of filter µ
##
plot_L_temp_m <- plotTemporal(odeLogistic, within(p_log, {r <- 0.4; trans <- 0.001; m <- 0.001}), filterslope = c(0.001, 0), times = t_temp,
                              title = "X) Temporally shifting realized space of carrying capacity niche.", plotregression = T)
##
plot_M_temp_m <- plotTemporal(odeMortality,  within(p_mort, {r <- 2; trans <- 0.0005; M_a <- 0.01}),
                              filter_j = c(0.4, 0.15), filter_a = c(0.4, 0.15), filterslope = c(0.0002, 0), times = t_fit,
                              title = "Temporally shifting environmental drivers;\n ontogenetically synchronous niche.",
                              plotregression = T, plotcenter = F, plotrange = F)

##
plot_M_both_m <- plotTemporal(odeMortality,  within(p_mort, {r <- 2; trans <- 0.0005; M_a <- 0.01}),
                              filter_j = c(0.4, 0.15), filter_a = c(0.49, 0.15), filterslope = c(0.0001, 0), times = t_fit,
                              title = "Temporally shifting environmental drivers\n and ontogenetic niche shift (opposing directions).",
                              plotregression = T, plotcenter = F, plotrange = F)



#### CASE Temp2): Temporal change of σ
## J < A
plot_L_temp_s <- plotTemporal(odeLogistic, within(p_log, {trans <- 0.001; m <- 0.001}), filterslope = c(0, -0.0002), times = t_temp,
                              title = "X) Temporally contracting realized space of carrying capacity niche", plotcenter = F)
# ## J > A works
# plot_L_temp_s <- plotTemporal(odeLogistic, within(p_log, {trans <- 0.001; m <- 0.001}), filterslope = c(0, 0.005), times = t_temp,
#                               title = "X) Temporally expanding realized space of carrying capacity niche", plotrange = T)

## J > A
plot_M_temp_s <- plotTemporal(odeMortality,  within(p_mort, {trans <- 0.001}), filterslope = c(0, 0.001), times = t_temp,
                              title = "X) Temporally expanding realized space of the mortality niche", plotcenter = F)

# ## J < A appears only achievable in edge cases with much higher juvenile mortality and historic effects
# plot_M_temp_s <- plotTemporal(odeMortality,  within(p_mort, {trans <- 0.02; M_j <- 0.4; M_a <- 0.1; d <- 0}), filterslope = c(0, -0.01), times = t_temp,
#                               title = "X) Temporally contracting realized space of the mortality niche", plotrange = T)


#### CASE: Trailing. Range expansion with trailing edge due to accumulating adults and pioneering juveniles.
## Trailing can be produced with both m and K
##

plot_L_temp_m_trailing <- plotTemporal(odeLogistic, within(p_log, {r = 0.4; m = 0.0002}), filterslope = c(0.001, 0), times = t_temp,
                              title = "Temporally shifting environmental driver;\n expansion at the trailing edge due to low mortality.", plottotal = T, plotcenter = F)

## 
plot_M_temp_m_trailing <- plotTemporal(odeMortality, within(p_mort, {r <- 1; trans <- 0.0028; M_a <- 0.005; M_j <- 0.1; d <- 50; K_j <- K_a <- 10000}),
                                       filterslope = c(0.001, 0), times = t_temp,
                                       filter_j = c(m = 0.4, s = 0.13),
                                       filter_a = c(m = 0.4, s = 0.13),
                                       title = "Temporally shifting environmental driver;\n expansion at the trailing edge due to low mortality.",
                                       plottotal = F, plotcenter = F)


## CASE: Ragging. Range contraction with dispersal limited leading edge.
## Expansion with trailing edge due to accumulating adults
# somepar <- within(p_mort, {trans <- 0.001; r <- 0.2; d <- 20}) ## Parameters for expansion because of slow growth

##
p_ragging <- list(r = 0.1, trans = 0.001, m = 0.006, K_j = 40, K_a = 10, d = 0, seedbank = 0, M_j = 0.01, M_a = 0.01)
plot_L_temp_m_ragging <- plotTemporal(odeLogistic, p_ragging, filterslope = c(0.002, 0), times = t_temp,
                                       title = "Temporally shifting environmental driver;\n", plottotal = T, plotcenter = F)

##
plot_M_temp_m_ragging <- plotTemporal(odeMortality,
                                       within(p_mort, {r <- 1; trans <- 0.0018; M_a <- 0.0001; M_j <- 0.6; d <- 1; K_j <- K_a <- 10000}),
                                       filterslope = c(0.001, 0), times = t_temp,
                                       filter_j = c(m = 0.4, s = 0.13),
                                       filter_a = c(m = 0.4, s = 0.13),
                                       title = "Temporally shifting environmental driver;\n contraction at the leading edge due to limited regeneration.",
                                      plottotal = F, plotcenter = F)

## Population growth lagging behind the filter shift can lead to two effects one after another:
## 1. Spreading occurrence, due to the the trailing edge in abundance still producing new abundance against increasing mortality.
## Only possible when abundance is limited by shifting K ("new opportunity").
## 2. Ragging occurrence because trailing dispersal is a bottleneck in providing K with new abundance at the leading filter edge.
## Improving dispersal will increase the spreading, but decrease the ragging effect.



####################################################################################
# Arrange plots of occurrence cases  ---------------------------------------------##
####################################################################################
blank <- ggplot() + theme_void()

legend <- simulateAbundance(odeMortality, p_mort, filterslope = c(0.000001, 0), times = c(1, 100))  %>%
  plotAbundance(legend = "b", plotfilter = "M", plottotal = T) %>%
  cowplot::get_legend()

legend_sum <- plotTemporal(odeLogistic, within(p_log, {r = 0.4; m = 0.0002}), filterslope = c(0.001, 0), times = t_temp,
                           plottotal = T, plotcenter = F, legend = "r") %>% 
  cowplot::get_legend()

## returns a ggplot object with the subheader removed and other labels replaced
rL <- function(ggplot, t = NULL, x = "environment at t = 1") ggplot + ggtitle(label = t, subtitle = NULL) + xlab(x)
rEnv <- function(ggplot, x = "environment at t = 1") ggplot + xlab(x)


## returns a text box with text c
ct <- function(c, ...) grobTree(rectGrob(), textGrob(c, ...))

## returns a list of the provided elements
itemize <- function(...) paste("- ", c(...), collapse = "\n")



# Plot grid 1: Regression -----------------------------------------------

shiftgrid <- cowplot::plot_grid(rEnv(plot_M_temp_m),
                                rEnv(plot_M_ont_m, x = "environment"),
                                rEnv(plot_M_both_m), legend_sum,
                                nrow = 2, ncol = 2,
                                labels = c("A", "B", "C", ""), label_size = 14)
cowplot::save_plot(glue("Publishing/Plots/Sim_shift.pdf"), shiftgrid, base_height = 9, base_asp = 1.45, device = cairo_pdf) # cairo_pdf will embed fonts!
gc()


# Plot grid 2. Effect of range on occurrence width and margins ----------------------------------------

plots_width <-  list(plot_M_same, legend_sum, plot_M_temp_m_trailing, plot_M_temp_m_ragging)
plotgrid_width <- cowplot::plot_grid(plotlist = plots_width, ncol = 2, nrow = 2,
                                     labels = c("D", "", "E", "F"), label_size = 14)
cowplot::save_plot(glue("Publishing/Plots/Sim_width.pdf"), plotgrid_width, base_height = 9, base_asp = 1.45, device = cairo_pdf) # cairo_pdf will embed fonts!



# Plot grid 1 and M_same combined -----------------------------------------------
shiftgrid_combined <- cowplot::plot_grid(rEnv(plot_M_temp_m),
                                         rEnv(plot_M_ont_m, x = "environment"),
                                         rEnv(plot_M_both_m),
                                         rEnv(plot_M_same), # get legend from somewhere else!
                                         nrow = 2, ncol = 2,
                                         labels = c("A", "B", "C", "D", ""), label_size = 14)
cowplot::save_plot(glue("Publishing/Plots/Sim.pdf"), shiftgrid_combined, base_height = 9, base_asp = 1.45, device = cairo_pdf) # cairo_pdf will embed fonts!
gc()



# # 2x2 Plotgrid 1. Ontogeny, time can have the same effects --------------------------------
# # occurrence cases:
# 
# shift_temp <- "shifting environmental driver" %>% ct()
# shift_ont <- "shifting niche" %>% ct()
# 
# JgA_temp <- "expanding environmental driver" %>% ct()
# JgA_ont <- "equal or contracting niche" %>% ct()
# 
# 
# 
# temprow <- grid.arrange(
#   ## Rows of the plot grid.
#   shift_temp, JgA_temp,
#   rL(plot_M_temp_m), rL(plot_M_temp_s),
# 
#   ### Layout matrix, arranged as in rows.
#   layout_matrix = rbind(c(1, 2),
#                         c(3, 4)),
#   heights = c(1, 3)
# )
# 
# ontrow <- grid.arrange(
#   shift_ont, JgA_ont,
#   rL(plot_M_ont_m, x = "environment"), rL(plot_M_same, x = "environment"),
#   
#   ### Layout matrix, arranged as in rows.
#   layout_matrix = rbind(c(1, 2),
#                         c(3, 4)),
#   heights = c(1, 3)
# )
# 
# gc()
# textframegrid_occ <- grid.arrange(blank, ct("Different occurence centers can be caused by …"), ct("Different occurrence widths can be caused by …"),
#                                   ct("temporal", rot = 60),                           temprow,                                        
#                                   ct("ontogenetic", rot = 60),                        ontrow,                                         
#                                  legend, textGrob("dotted lines are niches"),
#                                  
#                                  layout_matrix = rbind(c(1, 2, 3),
#                                                        c(4, 5, 5),
#                                                        c(6, 7, 7),
#                                                        c(8, 8, 9)),
#                                  heights = c(1, 4, 4, 0.5), widths = c(1, 3, 3))



# # 2x3 Plotgrid. Ontogeny, time can have the same effects --------------------------------
# 
# # occurrence cases:
# 
# shift_temp <- itemize("shifting realized space of K niche", "shifting realized space of m niche") %>% ct()
# shift_ont <- itemize("shifting K niche", "shifting m niche") %>% ct()
# 
# JsA_temp <- itemize("contracting realized space of K niche", "(fringe cases with m niche)") %>% ct()
# JsA_ont <- itemize("expanding K niche", "(fringe cases with m niche)") %>% ct() # extreme cases of much higher juvenile mortality with m niche
# 
# JgA_temp <- itemize("expanding realized space of K niche", "expanding realized space of m niche") %>% ct()
# JgA_ont <- itemize("contracting K niche", "contracting m niche", "equal m niche") %>% ct()
# 
# 
# 
# temprow <- grid.arrange(
#   ## Rows of the plot grid.
#   shift_temp, JsA_temp, JgA_temp,
#   plot_L_temp_m, plot_L_temp_s, plot_M_temp_s,
#   
#   ### Layout matrix, arranged as in rows.
#   layout_matrix = rbind(c(1, 2, 3),
#                         c(4, 5, 6)),
#   heights = c(1, 3)
# )
# 
# ontrow <- grid.arrange(
#   shift_ont, JsA_ont, JgA_ont,
#   plot_M_ont_m, plot_L_ont_s, plot_M_same, # 4, 5, 6
#   
#   ### Layout matrix, arranged as in rows.
#   layout_matrix = rbind(c(1, 2, 3),
#                         c(4, 5, 6)),
#   heights = c(1, 3)
# )
# 
# textframegrid_occ <- grid.arrange(ct("occurence shift"), ct("occurrence J < A"), ct("occurrence J > A"), blank,
#                                   temprow,                                        ct("temporal", rot = 310),
#                                   ontrow,                                         ct("ontogenetic", rot = 310),
#                                   legend, textGrob("dotted lines are niches,\n arrows correspond to 95% iqr"),
#                                   
#                                   layout_matrix = rbind(c(1, 2, 3, 4),
#                                                         c(5, 5, 5, 6),
#                                                         c(7, 7, 7, 8),
#                                                         c(9, 9, 10, 4)),
#                                   heights = c(1, 4, 4, 0.5))





####################################################################################
# Parameters table                   ---------------------------------------------##
####################################################################################

plots <- list(plot_M_temp_m, plot_M_ont_m, plot_M_both_m, plot_M_same)
plotnames <- c('plot_M_temp_m', "plot_M_ont_m", 'plot_M_both_m', 'plot_M_same')

## Former version with more plots
# plots <- list(plot_M_ont_m, plot_M_temp_m, plot_M_both_m, plot_M_same, plot_M_temp_m_trailing, plot_M_temp_m_ragging)
# plotnames <- c("plot_M_ont_m", 'plot_M_temp_m', 'plot_M_both_m', 'plot_M_same', 'plot_M_temp_m_trailing', 'plot_M_temp_m_ragging')

extractParameters <- function(plot) {
  parametername <- union(names(p_log), names(p_mort))
  par <- unlist(plot$simparameters)[parametername]
  filterpar <- plot$filterparameters
  time <-  unique(plot$A$time)
  if(length(time) == 1) time <- c(time1 = time, time2 = NA)
  c(par, filterpar, t = time)
}

Parameters <- t(sapply(plots, extractParameters))
rownames(Parameters) <- c("A temporal shift", "B ontongenetic niche shift", "C temp. and ont. niche shift", "D identical ont. niches") # dput(plotnames)

# xtable::xtable(Parameters)


opencurly <- "{"
closecurly <- "}"

Parameters <- as.data.frame(Parameters)

Partable_latex <- Parameters %>%
  as.data.frame() %>%

  ## This is just for reordering (Make sure that parameters are complete!):
  dplyr::select("$r$" = r, "$g$" = trans, "$K$" = K_j, "$\\hat{K}$" = K_a, "$M$" = M_j, "$\\hat{M}$" = M_a, "$d$"= d,
                "$a_{\\mu}$" = m_j, "$\\sigma$" = s_j, "$\\hat{a}_{\\mu}$" = m_a, "$\\hat{\\sigma}$" = s_a, "$b_{\\mu}$" = b_m,
                "$t_1$" = t1, "$t_2$" = t2) %>% # b_s
  kable("latex", caption = glue("Simulation parameters. \\label{opencurly}tab:parameters{closecurly}"),
        escape = F, booktabs = T) %>%
  add_header_above(c(" " = 3, "carrying cap." = 2, "max. mort." = 2, " " = 1, "f properties" = 5, "times" = 2)) %>%
  kable_styling(font_size = 7)

# cat(Partable_latex, file = glue("~/Documents/Studium/Projects/Papers/Mismatch/Tables/Partable_new.tex"))


