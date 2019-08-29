#' ---
#' title: "Bat Viral Detection and Reproduction Analysis - Plotting"
#' author: "EcoHealth Alliance M&A Team - Code Drafted by Evan Eskew"
#' ---

# /*
#==============================================================================
# */


#+ loading_chunk, echo=FALSE, message=FALSE, results="hide"


library(tidyverse)
library(rethinking)
library(assertthat)
library(ggridges)
library(viridis)
library(cowplot)
library(diagram)
library(maptools)
library(reskew) # devtools::install_github("eveskew/reskew")

logistic <- rethinking::logistic

# Load data

d.sample.sizes <- read.csv("outputs/d.sample.sizes.csv")
dat.f.trim <- readRDS("data/cleaned_data/dat.f.trim.rds")
dat.f.stan <- readRDS("stan/cleaned_data/dat.f.stan.rds")

# Load and process Stan models

pars.to.trim <- c(
  "alpha_host_species_offset_tilde",
  "alpha_year_tilde",
  "alpha_country_tilde",
  "alpha_specimen_type_group_tilde",
  "alpha_test_requested_mod_tilde",
  "alpha_diagnostic_laboratory_name_tilde",
  "alpha",
  "log_lik",
  "lp__"
)

model.names <- str_replace(list.files("stan/saved_models/"), ".rds", "")

for (model.name in model.names) {
    
    processed.model.name <- paste0(model.name, ".p")
    
    assign(model.name, 
           readRDS(paste0("stan/saved_models/", model.name, ".rds"))
    )
    
    assign(processed.model.name, 
           process_stanfit(get(model.name), pars.to.trim = pars.to.trim))
    
    print(paste("Processed model:", model.name))
    print(paste("Divergences:", toString(get(processed.model.name)$divergences)))
    print("Rhat summary:")
    print(get(processed.model.name)$Rhat.summary)
}

# Load labelling info

var.effect.group.labels <- c(
  "host_species_offset" = "Host Species",
  "year" = "Year of Sample Collection",
  "country" = "Country of Sample Collection", 
  "specimen_type_group" = "Specimen Type", 
  "test_requested_mod" = "Viral Test Protocol",
  "diagnostic_laboratory_name" = "Diagnostic Laboratory Conducting Testing"
)

# Set plotting defaults

custom_theme <- theme_minimal() + 
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    text = element_text(size = 14, color = "black")
  )

# /*
#==============================================================================
# */


# Figure 1 - Observed data summary


species.to.keep <- dat.f.trim %>%
  distinct(binomial, pregnant_mod) %>%
  group_by(binomial) %>%
  summarize(n = n()) %>%
  filter(n == 2) %>%
  pull(binomial)

helper.data <- dat.f.trim %>%
  filter(binomial %in% species.to.keep) %>%
  distinct(binomial, pregnant_mod) %>%
  mutate(color = "black") %>%
  arrange(binomial)

plot1 <- dat.f.trim %>%
  filter(binomial %in% species.to.keep) %>%
  ggplot(aes(x = pregnant_mod, y = virus_detected)) +
  coord_cartesian(xlim = c(-0.2, 1.2), ylim = c(0, 0.20)) +
  xlab("") + ylab("Observed Viral Prevalance") +
  scale_x_continuous(
    labels = c("Not Pregnant", "Pregnant"), 
    breaks = c(0, 1)
  ) +
  stat_summary(aes(group = binomial), 
               fun.y = "mean", 
               geom = "point",
               color = helper.data$color) +
  stat_summary(aes(group = binomial), 
               fun.y = "mean", 
               geom = "line", size = 1,
               color = alpha(helper.data$color, 0.3)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    text = element_text(size = 20, color = "black"),
    strip.text = element_text(face = "bold")
  )

dat.f.trim %>%
  filter(binomial %in% species.to.keep) %>%
  group_by(binomial, pregnant_mod) %>%
  summarize(
    n = n(),
    n_positive = sum(virus_detected),
    prev = n_positive/n*100
  ) %>%
  select(binomial, pregnant_mod, prev) %>%
  spread(key = pregnant_mod, value = prev) %>%
  mutate(pregnant_prev_lower = `1` <= `0`)
    
species.to.keep <- dat.f.trim %>%
  distinct(binomial, lactating_mod) %>%
  group_by(binomial) %>%
  summarize(n = n()) %>%
  filter(n == 2) %>%
  pull(binomial)

helper.data <- dat.f.trim %>%
  filter(binomial %in% species.to.keep) %>%
  distinct(binomial, lactating_mod) %>%
  mutate(color = "black") %>%
  arrange(binomial)

plot2 <- dat.f.trim %>%
  filter(binomial %in% species.to.keep) %>%
  ggplot(aes(x = lactating_mod, y = virus_detected)) +
  coord_cartesian(xlim = c(-0.2, 1.2), ylim = c(0, 0.20)) +
  xlab("") + ylab("Observed Viral Prevalance") +
  scale_x_continuous(
    labels = c("Not Lactating", "Lactating"), 
    breaks = c(0, 1)
  ) +
  stat_summary(aes(group = binomial), 
               fun.y = "mean", 
               geom = "point",
               color = helper.data$color) +
  stat_summary(aes(group = binomial), 
               fun.y = "mean", 
               geom = "line", size = 1,
               color = alpha(helper.data$color, 0.3)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    text = element_text(size = 20, color = "black"),
    strip.text = element_text(face = "bold")
  )

dat.f.trim %>%
  filter(binomial %in% species.to.keep) %>%
  group_by(binomial, lactating_mod) %>%
  summarize(
    n = n(),
    n_positive = sum(virus_detected),
    prev = n_positive/n*100
  ) %>%
  select(binomial, lactating_mod, prev) %>%
  spread(key = lactating_mod, value = prev) %>%
  mutate(lactating_prev_lower = `1` <= `0`)


plot_grid(plot1, plot2, nrow = 1, scale = 0.95,
          labels = "auto", label_size = 24)

ggsave("outputs/Fig1.png", width = 10, height = 6, dpi = 350)

# /*
#==============================================================================
# */


# Figure 2 - Model parameter summary


model.cols <- c("black", "deepskyblue3", "darkseagreen4", 
                "coral3", "plum4", "darkorange3")

pars <- c("mu_alpha", "beta_pregnant_mod", "beta_lactating_mod")
bracket <- 
  list(c("Reproductive Effects on Viral Detection", 
         "Pregnancy Effect", 
         "Lactation Effect")
  )

estimate.method <- "median"
conf.int <- TRUE
conf.method <- "HPDinterval"
conf.level <- 0.95


tidy.model.output <- data.frame(NULL)

for (model.name in model.names) {
  
  temp.df <- 
    broom::tidy(get(model.name), 
                estimate.method = estimate.method, conf.int = conf.int,
                conf.method = conf.method, conf.level = conf.level,
                pars = pars) %>%
    mutate(model = rep(model.name, length(pars)),
           model = str_replace(model, "model.f", ""),
           model = str_replace(model, ".", ""),
           model = ifelse(model == "", "All Viral Families", model))
  
  tidy.model.output <- rbind(tidy.model.output, temp.df)
}

helper.data <- dat.f.trim %>%
  group_by(test_requested_viral_family) %>%
  count()

helper.data <- bind_rows(
  helper.data, 
  data.frame(test_requested_viral_family = "All Viral Families", n = nrow(dat.f.trim)
  )
)

tidy.model.output <- tidy.model.output %>%
  left_join(., helper.data, 
            by = c("model" = "test_requested_viral_family")
  ) %>%
  mutate(model = paste0(model, " (n = ", n, ")"))

models.excluding.all <- 
  unique(tidy.model.output$model)[!str_detect(unique(tidy.model.output$model), "All Viral Families")]
all.model <- 
  unique(tidy.model.output$model)[str_detect(unique(tidy.model.output$model), "All Viral Families")]
rev.models.excluding.all <- rev(sort(models.excluding.all))

tidy.model.output <- tidy.model.output %>%
  mutate(
    model = factor(model, levels = c(rev.models.excluding.all, all.model))
  )


plot1 <- {
  tidy.model.output %>%
    dotwhisker::relabel_predictors(
      c(mu_alpha = "Global Intercept",
        beta_pregnant_mod = "Pregnancy Effect",
        beta_lactating_mod = "Lactation Effect"
      )
    ) %>%
    dotwhisker::dw_plot(
      dot_args = list(size = 2),
    ) +
    custom_theme +
    xlab("Parameter Estimate") +
    theme(legend.text = element_text(size = 8),
          legend.title = element_text(size = 9)) +
    geom_vline(xintercept = 0, colour = "black", linetype = 2) +
    xlim(-13, 5) +
    scale_color_manual(values = rev(model.cols),
                       name = "Viral Dataset") +
    guides(color = guide_legend(reverse = TRUE))
  } %>%
  dotwhisker::add_brackets(bracket)

plot2 <- model.f %>% 
  broom::tidy(., 
              estimate.method = estimate.method, conf.int = conf.int,
              conf.method = conf.method, conf.level = conf.level,
              pars = c("alpha_host_species")) %>%
  dotwhisker::dw_plot(.,
                      dot_args = list(col = "black"),
                      whisker_args = list(col = alpha("black", 0.8))
  ) +
  custom_theme +
  theme(axis.text.y = element_blank()) +
  xlab("Parameter Estimate") +
  geom_vline(xintercept = 0, colour = "black", linetype = 2) +
  xlim(-10, 1) +
  ggtitle("Intercepts by Host Species") +
  theme(plot.title = element_text(hjust = 0.5))

model.f %>% 
  broom::tidy(., 
              estimate.method = estimate.method, conf.int = conf.int,
              conf.method = conf.method, conf.level = conf.level,
              pars = c("alpha_host_species")) %>%
  mutate(implied_probability = logistic(estimate)) %>%
  arrange(implied_probability)


plot_grid(plot1, plot2, ncol = 1, scale = c(1, 0.95), rel_heights = c(1, 1.1),
          labels = "auto", label_size = 24)

ggsave("outputs/Fig2.png", width = 8, height = 10, dpi = 350)

plot1

ggsave("outputs/Fig2.png", width = 8, height = 6, dpi = 350)
 
# /*
#==============================================================================
# */


# Figure 3 - Model-based predictions


sims <- 50
samples.per.sim <- 1000


big.sim.df <- data.frame(NULL)

for(x in model.names) {
  
  dat <- get(paste0(x, ".p"))$df
  
  sim.df <- data.frame(
    rep(x, sims*3),
    rep(c("Non-reproductive", "Pregnant", "Lactating"), sims),
    rep(NA, sims*3),
    rep(NA, sims*3)
  )
  colnames(sim.df) <- c("model", "condition", "positives", "prevalence")
  sim.df
  
  set.seed(8)
  
  for (i in (1:sims)*3 - 2) {
    
    sim.df$positives[i] <-
      sum(rbinom(samples.per.sim, 
                 prob = sample(logistic(dat$mu_alpha), samples.per.sim), 
                 size = 1))
    
    sim.df$positives[i + 1] <-
      sum(rbinom(samples.per.sim, 
                 prob = sample(logistic(dat$mu_alpha + dat$beta_pregnant_mod), 
                               samples.per.sim), 
                 size = 1))
    
    sim.df$positives[i + 2] <-
      sum(rbinom(samples.per.sim, 
                 prob = sample(logistic(dat$mu_alpha + dat$beta_lactating_mod),
                               samples.per.sim), 
                 size = 1))
    
    sim.df$prevalence[i:(i + 2)] <- sim.df$positives[i:(i + 2)]/samples.per.sim*100
  }
  
  big.sim.df <- bind_rows(big.sim.df, sim.df)
}

big.sim.df$condition <- 
  factor(big.sim.df$condition, 
         levels = c("Non-reproductive", "Pregnant", "Lactating"))

model.labels <- list(
  'model.f' = "All Viral Families",
  'model.f.Adenoviridae' = "Adenoviridae",
  'model.f.Coronaviridae' = "Coronaviridae",
  'model.f.Herpesviridae' = "Herpesviridae",
  'model.f.Paramyxoviridae' = "Paramyxoviridae",
  'model.f.Polyomaviridae' = "Polyomaviridae"
)

model_labeller <- function(variable, value) {return(model.labels[value])}

alpha <- 1
color.values <- c(
  alpha("cornsilk3", alpha),
  alpha("darkseagreen", alpha),
  alpha("lightsteelblue", alpha)
)

plot <- big.sim.df %>%
  ggplot(aes(x = condition, y = prevalence, color = condition)) +
  geom_jitter(height = 0) +
  ylab("Predicted Viral Prevalence") +
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        text = element_text(size = 18, color = "black"),
        strip.text = element_text(size = 10, face = "bold"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal"
  ) +
  scale_color_manual(
    values = color.values, 
    name = "Reproductive Condition",
    aesthetics = c("colour", "fill")
  ) +
  facet_wrap(~model, labeller = model_labeller, nrow = 1)
  
summary.data <- big.sim.df %>%
  group_by(model, condition) %>%
  summarize(prevalence = mean(prevalence)) %>%
  mutate(x = rep(1:3)) %>%
  ungroup()


plot + 
  geom_path(
    data = summary.data, inherit.aes = FALSE,
    aes(x = x, y = prevalence),
    color = alpha("black", 0.8), 
    size = 0.7,
    arrow = arrow(angle = 20, length = unit(0.1, "inches"), type = "closed")
  )

ggsave("outputs/Fig3.png", height = 5, width = 8, dpi = 350)

# /*
#==============================================================================
# */


# Supplementary Table 1


dat.f.trim %>%
  group_by(binomial) %>%
  summarize(
    "Viral Tests" = n(),
    "Positive Viral Tests" = sum(virus_detected == 1),
    "Viral Tests from Pregnant Individuals" = sum(pregnant_mod == 1),
    "Viral Tests from Lactating Individuals" = sum(lactating_mod == 1)
  ) %>%
  mutate(binomial = as.character(binomial)) %>%
  dplyr::rename("Host Species" = binomial) %>%
  arrange(`Host Species`) %>%
  ungroup() %>%
  rbind(
    c("Total",
      sum(.$"Viral Tests", na.rm = TRUE),
      sum(.$"Positive Viral Tests", na.rm = TRUE),
      sum(.$"Viral Tests from Pregnant Individuals", na.rm = TRUE),
      sum(.$"Viral Tests from Lactating Individuals", na.rm = TRUE)
    )
  ) %>%
  mutate_all(funs(replace(., is.na(.), " "))) %>%
  write_csv(., "outputs/TableS1.csv")

# /*
#==============================================================================
# */


# Figure S1 - Data filtering and sample sizes flowchart


my.labels <- c(
  paste0("All PREDICT-1 Data (in EIDITH Database)\nNumber of individual animals: ", 
         d.sample.sizes[which(d.sample.sizes == "d2"), 2],
         ", Number of species: ",
         d.sample.sizes[which(d.sample.sizes == "d2"), 3]), 
  'Filter for:\norder == "Chiroptera"',
  paste0("PREDICT-1 Bat Dataset\nNumber of individual animals: ",
         d.sample.sizes[which(d.sample.sizes == "d.bat"), 2],
         ", Number of species: ",
         d.sample.sizes[which(d.sample.sizes == "d.bat"), 3]),
  'Filter for:\npredict_protocol == "TRUE", animal_classification == "Wild",
  age_class is known, non-pooled cPCR tests, testing results available',
  paste0("Quality-filtered PREDICT-1 Bat Dataset\nNumber of individual animals: ",
         d.sample.sizes[which(d.sample.sizes == "d3"), 2],
         ", Number of species: ",
         d.sample.sizes[which(d.sample.sizes == "d3"), 3]),
  'Filter for:\nage_class %in% c("Subadult", "Adult"), female animals, 
  reproductive data available, exclusion of viral groups without positives',
  paste0("Final Female Bat Dataset\nNumber of individual animals: ",
         d.sample.sizes[which(d.sample.sizes == "dat.f.trim"), 2],
         ", Number of species: ",
         d.sample.sizes[which(d.sample.sizes == "dat.f.trim"), 3])
)

my.text.size <- 1
rx <- 0.35
ry <- 0.05
lwd <- 1
arr.pos <- 0.48
box.col1 <- "rosybrown1"
box.col2 <- "palegreen2"
round.col <- "snow1"
shadow.col <- alpha("white", 0.0)


png("outputs/FigS1.png", width = 1400, height = 1200, res = 160)

par(mar = c(0, 0, 0, 0) + 0.05,
    family = "Arial")

openplotmat()

pos <- diagram::coordinates(c(1, 1, 1, 1, 1, 1, 1))
straightarrow(from = pos[1, ], to = pos[2, ], lwd = lwd, arr.pos = arr.pos)
straightarrow(from = pos[2, ], to = pos[3, ], lwd = lwd, arr.pos = arr.pos)
straightarrow(from = pos[3, ], to = pos[4, ], lwd = lwd, arr.pos = arr.pos)
straightarrow(from = pos[4, ], to = pos[5, ], lwd = lwd, arr.pos = arr.pos)
straightarrow(from = pos[5, ], to = pos[6, ], lwd = lwd, arr.pos = arr.pos)
straightarrow(from = pos[6, ], to = pos[7, ], lwd = lwd, arr.pos = arr.pos)

textrect(mid = pos[1, ], lab = my.labels[1], cex = my.text.size, 
         radx = rx, rady = ry, box.col = box.col1, shadow.col = shadow.col)
textround(mid = pos[2, ], lab = my.labels[2], cex = my.text.size, 
          radx = rx, rady = ry, box.col = round.col, shadow.col = shadow.col)
textrect(mid = pos[3, ], lab = my.labels[3], cex = my.text.size, 
         radx = rx, rady = ry, box.col = box.col1, shadow.col = shadow.col)
textround(mid = pos[4, ], lab = my.labels[4], cex = my.text.size, 
          radx = rx, rady = ry, box.col = round.col, shadow.col = shadow.col)
textrect(mid = pos[5, ], lab = my.labels[5], cex = my.text.size, 
         radx = rx, rady = ry, box.col = box.col1, shadow.col = shadow.col)
textround(mid = pos[6, ], lab = my.labels[6], cex = my.text.size, 
          radx = rx, rady = ry, box.col = round.col, shadow.col = shadow.col)
textrect(mid = pos[7, ], lab = my.labels[7], cex = my.text.size, 
         radx = rx, rady = ry, box.col = box.col2, shadow.col = shadow.col)

dev.off()

# /*
#==============================================================================
# */


# Figure S3 - Sampling map


summary.mapping.df <- dat.f.trim %>%
  group_by(country) %>%
  summarize(n = n_distinct(test_id)) %>%
  mutate(
    country = case_when(
      country == "Tanzania" ~ "United Republic of Tanzania",
      country == "Vietnam" ~ "Viet Nam",
      TRUE ~ country
    )
  )

color.df <- 
  data_frame(
    country = wrld_simpl@data$NAME,
    color = rep(alpha("darkgreen", 0.3), length(wrld_simpl@data$NAME))
  ) %>%
  left_join(., summary.mapping.df, by = "country") %>%
  mutate(
    color = ifelse(
      !is.na(n), alpha("firebrick", logistic(log(n/100))), color
    )
  )

legend.df <-
  data_frame(
    n_tests = c(1, 10, 100, 1000, 10000)
  ) %>%
  mutate(
    fill = alpha("firebrick", logistic(log(n_tests/100)))
  )

data("wrld_simpl")


png(filename = "outputs/FigS3.png", width = 1000, height = 600, res = 80)

plot(wrld_simpl, xlim = c(-120, 150), ylim = c(-60, 60),
     col = color.df$color,
     bg = alpha("lightblue1", 0.05))

legend(x = -125, y = -10, 
       legend = legend.df$n_tests, fill = legend.df$fill,
       title = "No. cPCR tests")

dev.off()

# /*
#==============================================================================
# */


# Figure S4 - Trace plots


p <- rstan::traceplot(
  model.f,
  pars = c("mu_alpha", 
           "beta_pregnant_mod", "beta_lactating_mod", 
           "sigma_vector", "scale_for_sigmas"),
  ncol = 3
)

levels(p$data$parameter) <- c(
  "Global Intercept", 
  "Pregnancy Effect", "Lactation Effect",
  "σ (Host Species)", "σ (Year)", "σ (Country)", 
  "σ (Specimen Type)", "σ (Viral Test Protocol)", "σ (Diagnostic Laboratory)",
  "σ (Hierarchical Scale)"
)


p + scale_color_viridis(discrete = TRUE, option = "plasma", end = 0.85) +
  ggtitle("Parameter trace plots for Bayesian model of viral detection in adult female bats") +
  theme(
    text = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 18),
    strip.text.x = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

ggsave("outputs/FigS4.png", width = 10, height = 12, dpi = 350)

# /*
#==============================================================================
# */


# Figure S5 - Ridgeline plots for varying effects


dat.df <- dat.f.trim
model.df <- model.f.p$df

cols.to.plot <- grep("alpha", colnames(model.df), value = T) %>%
  grep("mu_|tilde|species\\[", ., value = T, invert = T)

var.effect.labels <- 
  list("host_species_offset" = sort(unique(dat.df$binomial)),
       "year" = levels(factor(dat.df$year)),
       "country" = levels(factor(dat.df$country)),
       "specimen_type_group" = sort(unique(dat.df$specimen_type_group)),
       "test_requested_mod" = sort(unique(dat.df$test_requested_mod)),
       "diagnostic_laboratory_name" = sort(unique(dat.df$diagnostic_laboratory_name))
  )

assert_that(length(cols.to.plot) == length(flatten(var.effect.labels)))


# To plot all varying effects

png("outputs/FigS5a.png", width = 1200, height = 800)

p <- model.df %>% 
  select(one_of(cols.to.plot)) %>%
  gather_("parameter", "value", cols.to.plot, factor_key = T) %>%
  mutate(
    varying_effect_group = gsub("\\[[0-9]+\\]", "", parameter) %>% gsub("alpha_", "", .),
    varying_effect_group = factor(varying_effect_group, names(var.effect.group.labels))
  ) %>%
  # To order the distributions by median values
  mutate(parameter = reorder(parameter, value, median)) %>%
  
  ggplot(aes(x = value, y = parameter, 
             height = ..density.., fill = varying_effect_group)) +
  xlab("Parameter Value") + ylab("Varying Intercepts Cluster") +
  geom_density_ridges() +
  xlim(-10, 10) +
  # To facet by varying effects group
  facet_wrap(~varying_effect_group, scale = "free_y", 
             labeller = as_labeller(var.effect.group.labels)) +
  theme_ridges(center_axis_labels = TRUE, font_size = 16) +
  # To remove y-axis labels
  theme(axis.text.x = element_text(size = 20),
        axis.text.y = element_blank(), 
        axis.title = element_text(size = 32),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 18, face = "bold", margin = margin(b = 5))) +
  coord_cartesian(clip = "off")

p

dev.off()


# To plot varying effects one at a time

plotting.list <- list(
  c("host_species_offset", "year", "country",
    "specimen_type_group", "test_requested_mod", "diagnostic_laboratory_name"),
  c("outputs/FigS5b.png", "outputs/FigS5c.png", "outputs/FigS5d.png", 
    "outputs/FigS5e.png", "outputs/FigS5f.png", "outputs/FigS5g.png"),
  c("Host Species", "Year of Sample Collection",
    "Country of Sample Collection", "Specimen Type", 
    "Viral Test Protocol", "Diagnostic Laboratory Conducting Testing"
     ),
  c(-6, -6, -6, -6, -10, -6)
)

for (i in 1:length(plotting.list[[1]])) {
  
  cols.to.plot <- grep(paste0(plotting.list[[1]][i], "\\["), 
                       colnames(model.df), value = T)
  
  labels <- var.effect.labels[[plotting.list[[1]][i]]]
  
  assert_that(length(cols.to.plot) == length(labels))
  
  fill.color <- ggplot_build(p)$data[[1]] %>%
    distinct(fill) %>%
    slice(i) %>%
    unlist()
  
  png(plotting.list[[2]][i], width = 1000, height = 1200)
  
  print(
    model.df %>% 
      select(one_of(cols.to.plot)) %>%
      gather_("parameter", "value", cols.to.plot, factor_key = T) %>%
      mutate(
        varying_effect_group = gsub("\\[[0-9]+\\]", "", parameter) %>% gsub("alpha_", "", .)
      ) %>%
      mutate(parameter = plyr::mapvalues(parameter, cols.to.plot, labels)) %>%
      # To order the distributions by median values
      mutate(parameter = reorder(parameter, value, median)) %>%
      
      ggplot(aes(x = value, y = parameter, 
                 height = ..density..)) +
      xlab("Parameter Value") + ylab(plotting.list[[3]][i]) +
      geom_density_ridges(fill = fill.color) +
      coord_cartesian(xlim = c(plotting.list[[4]][i], -1*plotting.list[[4]][i])) +
      theme_ridges(center_axis_labels = TRUE, font_size = 12) +
      theme(axis.title = element_text(size = 32),
            axis.text.x = element_text(size = 20),
            axis.text.y = element_text(size = 16)) +
      scale_y_discrete(expand = expand_scale(add = c(0.2, 1.5)))
  )
  
  dev.off()
}
