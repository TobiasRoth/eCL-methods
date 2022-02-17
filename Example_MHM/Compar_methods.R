rm(list = ls(all = TRUE))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Settings and load data----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Libraries
library(tidyverse)
library(patchwork)
library(ggthemes)
library(arm)
library(rjags)
library(TITAN2)

# Load data
dat <- read_csv("Example_MHM/Site_data.csv")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# (a) und (e) Visual inspection of gradient categories  ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Split gradient into categories
dat$NKat <- "<7.5"
dat$NKat[dat$Ndep >= 7.5] <- "7.5-12.5"
dat$NKat[dat$Ndep >= 12.5] <- "12.5-17.5"
dat$NKat[dat$Ndep >= 17.5] <- "17.5-22.5"
dat$NKat[dat$Ndep >= 22.5] <- ">22.5"
dat$NKat <- factor(dat$NKat, levels = c("<7.5", "7.5-12.5", "12.5-17.5", "17.5-22.5", ">22.5"))
table(dat$NKat)

# Make plot with total species richness
pd <- tibble(
  x = c(2, 3, 3, 2),
  y = c(0, 0, 60, 60)
)
plota <- dat %>% 
  group_by(NKat) %>% 
  dplyr::summarise(
    N = n(),
    meanSR = mean(SR),
    lo = t.test(SR)$conf.int[1],
    up = t.test(SR)$conf.int[2]
  ) %>% 
  ggplot(aes(x = NKat, y = meanSR)) +
  geom_point() +
  geom_errorbar(aes(ymin = lo, ymax = up), width = 0.2) +
  labs(
    title = "Visual inspection of gradient categories\n",
    subtitle = "(a) Estimated critical load: ~12.5 kg", 
    x = "",
    y = "Total number of species"
  ) +
  ylim(0, 60) +
  geom_polygon(data = pd, mapping = aes(x = x, y=y), alpha = 0.2, fill = "orange") 

# Make plot with oligotrophic species
pd <- tibble(
  x = c(2, 3, 3, 2),
  y = c(0, 0, 30, 30)
)
plote <- dat %>% 
  group_by(NKat) %>% 
  dplyr::summarise(
    N = n(),
    meanSR = mean(SR_oligo),
    lo = t.test(SR_oligo)$conf.int[1],
    up = t.test(SR_oligo)$conf.int[2]
  ) %>% 
  ggplot(aes(x = NKat, y = meanSR)) +
  geom_point() +
  geom_errorbar(aes(ymin = lo, ymax = up), width = 0.2) +
  labs(
    subtitle = "(e) Estimated critical load: ~12.5 kg", 
    x = "",
    y = "Number of oligotrophic species"
  ) +
  ylim(0, 30) +
  geom_polygon(data = pd, mapping = aes(x = x, y=y), alpha = 0.2, fill = "orange") 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# (b) TITAN approach: total numnber species  ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Load and prepare site x species matrix
getpldat <- function(s) {
  res <- numeric(nrow(dat))
  res[pl$siteID[pl$aID_SP == s]]<- pl$Deckung[pl$aID_SP == s]
  res
}
pl <- read_csv("Example_MHM/plant_data.csv")
species <- unique(pl$aID_SP) %>% sort()
speciesmat <- map_dfc(species, getpldat)
names(speciesmat) <- paste0("S", species)
speciesmat <- speciesmat[, apply(speciesmat, 2, function(x) sum(x > 0)) > 5] 

# Run titan
# titan_all <-
#   titan(
#     txa = speciesmat,
#     env = dat$Ndep,
#     numPerm = 100, # etwa 500 für die echte Analyse
#     nBoot = 500, # 500 oder 1000 für die echte Analyse
#     ncpus = 4)
# save(titan_all, file = "Example_MHM/titan_all.RData")
load("Example_MHM/titan_all.RData")
titan_all$sumz.cp["sumz-", c("cp", "0.05", "0.95")] %>% round(1)

# Make plot
d <- titan_all[["sppmax"]] %>% 
  as_tibble() %>% 
  mutate(aID_SP = names(speciesmat)) %>% 
  filter(filter == 1) %>% 
  mutate(aID_SP = fct_reorder(aID_SP, `50%`))
pd <- tibble(
  x = c(titan_all$sumz.cp["sumz-", "0.05"], titan_all$sumz.cp["sumz-", "0.95"], titan_all$sumz.cp["sumz-", "0.95"], titan_all$sumz.cp["sumz-", "0.05"]),
  y = c(0, 0, nrow(d) + 1, nrow(d) + 1)
)
plotb <- d %>% 
  ggplot(aes(y = aID_SP, x = `50%`)) +
  geom_point(col = "grey", cex = 0.3) +
  geom_errorbar(aes( xmin = `5%`, xmax = `90%`), width = 0, col = "grey") +
  geom_polygon(data = pd, mapping = aes(x = x, y=y), alpha = 0.2, fill = "orange") +
  geom_vline(xintercept = titan_all$sumz.cp["sumz-", "0.50"], col = "orange") +
  xlim(0, 40) +
  labs(
    title = "Threshold Indicator Taxa Analysis\n(TITAN)",
    subtitle = paste("(b) Estimated critical load: ", titan_all$sumz.cp["sumz-", "0.50"] %>% round(1), "kg"), 
    x = "",
    y = "Total number of species"
  ) +
  theme(axis.text.y = element_blank())

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# (f) TITAN approach: oligotrophic species  ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Load and prepare site x species matrix
getpldat <- function(s) {
  res <- numeric(nrow(dat))
  res[pl$siteID[pl$aID_SP == s]]<- pl$Deckung[pl$aID_SP == s]
  res
}
pl <- read_csv("Example_MHM/plant_data.csv") %>% 
  filter(N <= 2)
species <- unique(pl$aID_SP) %>% sort()
speciesmat <- map_dfc(species, getpldat)
names(speciesmat) <- paste0("S", species)
speciesmat <- speciesmat[, apply(speciesmat, 2, function(x) sum(x > 0)) > 5] 

# Run titan
# titan_oligo <-
#   titan(
#     txa = speciesmat,
#     env = dat$Ndep,
#     numPerm = 100, # etwa 500 für die echte Analyse
#     nBoot = 500, # 500 oder 1000 für die echte Analyse
#     ncpus = 4)
# save(titan_oligo, file = "Example_MHM/titan_oligo.RData")
load("Example_MHM/titan_oligo.RData")
titan_oligo$sumz.cp["sumz-", c("cp", "0.05", "0.95")] %>% round(1)

# Make plot
d <- titan_oligo[["sppmax"]] %>% 
  as_tibble() %>% 
  mutate(aID_SP = names(speciesmat)) %>% 
  filter(filter == 1) %>% 
  mutate(aID_SP = fct_reorder(aID_SP, `50%`))
pd <- tibble(
  x = c(titan_oligo$sumz.cp["sumz-", "0.05"], titan_oligo$sumz.cp["sumz-", "0.95"], titan_oligo$sumz.cp["sumz-", "0.95"], titan_oligo$sumz.cp["sumz-", "0.05"]),
  y = c(0, 0, nrow(d) + 1, nrow(d) + 1)
)
plotf <- d %>% 
  ggplot(aes(y = aID_SP, x = `50%`)) +
  geom_point(col = "grey") +
  geom_errorbar(aes( xmin = `5%`, xmax = `90%`), width = 0, col = "grey") +
  geom_polygon(data = pd, mapping = aes(x = x, y=y), alpha = 0.2, fill = "orange") +
  geom_vline(xintercept = titan_oligo$sumz.cp["sumz-", "0.50"], col = "orange") +
  xlim(0, 40) +
  labs(
    subtitle = paste("(f) Estimated critical load: ", sprintf( "%.1f", titan_oligo$sumz.cp["sumz-", "0.50"]), "kg"), 
    x = "",
    y = "Number of oligotrophic species"
  ) +
  theme(axis.text.y = element_blank())

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# (c) und (g) Point at which significance reduction can be observed  ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Prepare data 
d <- dat %>% 
  mutate(
    N = Ndep,
    ele = (elevation - 1500) / 200,
    incli = (inclination - 10) / 10,
    preci = (precipitation - 1500) / 1000,
    CACO3 = CACO3 - 3,
    hum = humidity - 3,
    light = dat$light - 3
  )

# Total number of species
mod <- glm(SR ~ Ndep + ele + incli + preci + CACO3 + hum, data = d, family = poisson) 
summary(mod)
nsim <- 1000
post <- sim(mod, nsim)@coef
lowintercept <-  post[, "(Intercept)"] %>% exp %>% quantile(probs = 0.05)
gradient <- seq(0, 40, 0.05) 
simres <- array(NA, dim = c(length(gradient), nsim))
for(i in 1:nsim) {
  simres[, i] <- exp(post[i, "(Intercept)"] + post[i, "Ndep"] * gradient)
}
predicted <- 
  tibble(
    med = apply(simres, 1, quantile, probs = 0.5),
    lo = apply(simres, 1, quantile, probs = 0.05),
    up = apply(simres, 1, quantile, probs = 0.95),
    NDep = gradient
  )
cl = gradient[which(apply(simres, 1, quantile, probs = 0.975) <= lowintercept) %>% min] 
plotc <- dat %>% 
  ggplot(aes(x = Ndep, y = SR)) +
  geom_point(cex = 0.6, col = "grey") +
  xlim(0, 40) +
  ylim(0, 80) +
  geom_line(data = predicted, aes(x = NDep, y = med), lwd = 0.8) +
  geom_line(data = predicted, aes(x = NDep, y = lo), lty = 2, lwd = 0.6) +
  geom_line(data = predicted, aes(x = NDep, y = up), lty = 2, lwd = 0.6) +
  geom_segment(aes(x = 0, y = lowintercept, xend = cl , yend = lowintercept), arrow = arrow(length = unit(0.3, "cm")), col = "orange") +
  geom_segment(aes(x = cl, y = lowintercept, xend = cl , yend = 0), arrow = arrow(length = unit(0.3, "cm")), col = "orange") +
  labs(
    title = "Point at which significance reduction\ncan be observed",
    subtitle = paste("(c) Estimated critical load:", cl, "kg"), 
    x = expression(paste("Nitrogen deposition [kg N ", ha^-1, " ", yr^-1, "]")),
    y = "Total number of species"
  ) 

# Oligotrophic species
mod <- glm(SR_oligo ~ Ndep + ele + incli + preci + CACO3 + hum, data = d, family = poisson) 
summary(mod)
nsim <- 1000
post <- sim(mod, nsim)@coef
lowintercept_oligo <-  post[, "(Intercept)"] %>% exp %>% quantile(probs = 0.05)
gradient <- seq(0, 40, 0.05) 
simres <- array(NA, dim = c(length(gradient), nsim))
for(i in 1:nsim) {
  simres[, i] <- exp(post[i, "(Intercept)"] + post[i, "Ndep"] * gradient)
}
predicted <- 
  tibble(
    med = apply(simres, 1, quantile, probs = 0.5),
    lo = apply(simres, 1, quantile, probs = 0.05),
    up = apply(simres, 1, quantile, probs = 0.95),
    NDep = gradient
  )
cl_oligo = gradient[which(apply(simres, 1, quantile, probs = 0.975) <= lowintercept_oligo) %>% min] 
plotg <- dat %>% 
  ggplot(aes(x = Ndep, y = SR_oligo)) +
  geom_point(cex = 0.6, col = "grey") +
  xlim(0, 40) +
  ylim(0, 50) +
  geom_line(data = predicted, aes(x = NDep, y = med), lwd = 0.8) +
  geom_line(data = predicted, aes(x = NDep, y = lo), lty = 2, lwd = 0.6) +
  geom_line(data = predicted, aes(x = NDep, y = up), lty = 2, lwd = 0.6) +
  geom_segment(aes(x = 0, y = lowintercept_oligo, xend = cl_oligo , yend = lowintercept_oligo), arrow = arrow(length = unit(0.3, "cm")), col = "orange") +
  geom_segment(aes(x = cl_oligo, y = lowintercept_oligo, xend = cl_oligo , yend = 0), arrow = arrow(length = unit(0.3, "cm")), col = "orange") +
  labs(
    subtitle = paste("(g) Estimated critical load:", cl_oligo %>% round(1), "kg"), 
    x = expression(paste("Nitrogen deposition [kg N ", ha^-1, " ", yr^-1, "]")),
    y = "Number of oligotrophic species"
  ) 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# (d) und (h) linear model with change-point  ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# JAGS Settings
t.n.thin <-      2
t.n.chains <-    2
t.n.burnin <- 5000
t.n.iter <-   2000

# Function to create initial values
inits <- function() {
  list(
    CL = rnorm(1, 15, 5),
    beta1 = rnorm(1, 0, 0.1),
    beta2 = rnorm(1, 0, 0.1),
    beta3 = rnorm(1, 0, 0.1),
    beta4 = rnorm(1, 0, 0.1),
    beta5 = rnorm(1, 0, 0.1),
    betaN = rnorm(1, 0, 0.1)
  )
}

# Prepare data for jags
nsite <- as.integer(nrow(dat))
N <- dat$Ndep
ele = (dat$elevation - 1500) / 200
incli = (dat$inclination - 10) / 10
preci = (dat$precipitation - 1500) / 1000
CACO3 = dat$CACO3 - 3
hum = (dat$humidity - 3)

# Make plot for total numbner of species
SR <- as.integer(dat$SR)
datjags <- list(SR = SR, nsite = nsite, N = N, ele = ele, incli = incli, preci = preci, CACO3 = CACO3, hum = hum)
jagres <- jags.model('Example_MHM/Change-point_model.R',data = datjags, n.chains = t.n.chains, inits = inits, n.adapt = t.n.burnin)
params <- c("CL", "SRpre")
mod <- coda.samples(jagres, params, n.iter=t.n.iter, thin=t.n.thin)
tymax <- 80
predicted <- summary(mod)$quantiles %>% 
  as_tibble(rownames = "variable") %>% 
  filter(str_detect(variable, "SRpre")) %>% 
  mutate(NDep = 1:40)
plotd <- dat %>% 
  ggplot(aes(x = Ndep, y = SR)) +
  geom_point(cex = 0.6, col = "grey") +
  xlim(0, 40) +
  geom_line(data = predicted, aes(x = NDep, y = `50%`), lwd = 0.8) +
  geom_line(data = predicted, aes(x = NDep, y = `2.5%`), lty = 2, lwd = 0.6) +
  geom_line(data = predicted, aes(x = NDep, y = `97.5%`), lty = 2, lwd = 0.6) +
  annotate(
    geom = "rect", 
    xmin = summary(mod)$quantiles["CL", "2.5%"], 
    xmax = summary(mod)$quantiles["CL", "97.5%"], 
    ymin = 0, ymax = tymax, alpha = .2,fill = "orange") +
  geom_line(
    data = tibble(
      x = rep(summary(mod)$quantiles["CL", "50%"], 2),
      y = c(0, tymax)), 
    aes(x = x, y = y),
    col = "orange") +
  labs(
    title = "Linear model with change-point\n",
    subtitle = paste("(d) Estimated critical load:", summary(mod)$quantiles["CL", "50%"] %>% round(1), "kg"), 
    x = expression(paste("Nitrogen deposition [kg N ", ha^-1, " ", yr^-1, "]")),
    y = "Total number of species"
  ) 

# Make plot for oligotrophic species
SR <- as.integer(dat$SR_oligo)
datjags <- list(SR = SR, nsite = nsite, N = N, ele = ele, incli = incli, preci = preci, CACO3 = CACO3, hum = hum)
jagres <- jags.model('Example_MHM/Change-point_model.R',data = datjags, n.chains = t.n.chains, inits = inits, n.adapt = t.n.burnin)
params <- c("CL", "SRpre")
mod <- coda.samples(jagres, params, n.iter=t.n.iter, thin=t.n.thin)
tymax <- 50
predicted <- summary(mod)$quantiles %>% 
  as_tibble(rownames = "variable") %>% 
  filter(str_detect(variable, "SRpre")) %>% 
  mutate(NDep = 1:40)
ploth <- dat %>% 
  ggplot(aes(x = Ndep, y = SR_oligo)) +
  geom_point(cex = 0.6, col = "grey") +
  xlim(0, 40) +
  geom_line(data = predicted, aes(x = NDep, y = `50%`), lwd = 0.8) +
  geom_line(data = predicted, aes(x = NDep, y = `2.5%`), lty = 2, lwd = 0.6) +
  geom_line(data = predicted, aes(x = NDep, y = `97.5%`), lty = 2, lwd = 0.6) +
  annotate(
    geom = "rect", 
    xmin = summary(mod)$quantiles["CL", "2.5%"], 
    xmax = summary(mod)$quantiles["CL", "97.5%"], 
    ymin = 0, ymax = tymax, alpha = .2,fill = "orange") +
  geom_line(
    data = tibble(
      x = rep(summary(mod)$quantiles["CL", "50%"], 2),
      y = c(0, tymax)), 
    aes(x = x, y = y),
    col = "orange") +
  labs(
    subtitle = paste("(h) Estimated critical load:", summary(mod)$quantiles["CL", "50%"] %>% round(1), "kg"), 
    x = expression(paste("Nitrogen deposition [kg N ", ha^-1, " ", yr^-1, "]")),
    y = "Number of oligotrophic species"
  ) 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Make figure ----
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# ggplot settings
theme_set(
  theme_clean() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      plot.subtitle = element_text(size = 10),
      plot.title.position = "plot",
      legend.title = element_blank(), 
      legend.position = "none", 
      legend.background = element_rect(colour = "white"),
      plot.background = element_blank())
)

# Add all plots to one figure
pdf("Example_MHM/Fig_MHM_method_comparison.pdf", width = 14, height = 8)
plota + plotb + plotc + plotd + plote + plotf + plotg + ploth +
  plot_layout(ncol = 4) &
  labs(x = expression(paste("Nitrogen deposition [kg N ", ha^-1, " ", yr^-1, "]")))
dev.off()
