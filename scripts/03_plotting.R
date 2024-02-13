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

data(wrld_simpl)
logistic <- rethinking::logistic

# Load data

d.sample.sizes <- read.csv("outputs/d.sample.sizes.csv")
dat.f.trim <- readRDS("data/cleaned_data/dat.f.trim.rds")
dat.f.stan <- readRDS("stan/cleaned_data/dat.f.stan.rds")

# Load and process Stan models

pars <- c("mu_alpha", "alpha_", "beta", "sigma")

model.names <- str_replace(list.files("stan/saved_models/"), ".rds", "")
processed.model.names <- paste0(model.names, ".p")

for (model.name in model.names) {
  
  # Load fit model
  assign(
    model.name, 
    readRDS(paste0("stan/saved_models/", model.name, ".rds"))
  )
  
  # Generate processed model with a subset of parameters in data frame format
  processed.model.name <- paste0(model.name, ".p")
  
  assign(
    processed.model.name,
    get(model.name)$draws() %>%
      posterior::as_draws_df() %>%
      select(contains(pars)) %>%
      select(!contains("tilde"))
  )
  
  print(paste("Processed model:", model.name))
  
  # Print per-chain divergent transitions using full fit model
  print(paste("Divergences:"))
  get(model.name)$diagnostic_summary("divergences", quiet = TRUE) %>%
    unlist() %>%
    as.vector() %>%
    print()
  
  # Print R-hat summary for all parameters in the processed model
  print("Rhat summary:")
  get(processed.model.name) %>%
    posterior::summarise_draws("rhat") %>%
    pull(rhat) %>%
    summary() %>%
    print()
}

# Load labeling info

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


# Which species have non-pregnant and pregnant individuals in the modeling 
# dataset?
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
  xlab("") + 
  ylab("Observed Viral Detection Probability") +
  scale_x_continuous(
    labels = c("Not Pregnant", "Pregnant"), 
    breaks = c(0, 1)
  ) +
  stat_summary(
    aes(group = binomial),
    fun = "mean",
    geom = "point",
    color = helper.data$color
  ) +
  stat_summary(
    aes(group = binomial),
    fun = "mean",
    geom = "line", linewidth = 1,
    color = alpha(helper.data$color, 0.3)
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    text = element_text(size = 20, color = "black"),
    strip.text = element_text(face = "bold")
  )

# Table summarizing the pregnancy dataset
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
    
# Which species have non-lactating and lactating individuals in the modeling 
# dataset?
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
  xlab("") + 
  ylab("Observed Viral Detection Probability") +
  scale_x_continuous(
    labels = c("Not Lactating", "Lactating"), 
    breaks = c(0, 1)
  ) +
  stat_summary(
    aes(group = binomial), 
    fun = "mean", 
    geom = "point",
    color = helper.data$color
  ) +
  stat_summary(
    aes(group = binomial), 
    fun = "mean", 
    geom = "line", linewidth = 1,
    color = alpha(helper.data$color, 0.3)
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    text = element_text(size = 20, color = "black"),
    strip.text = element_text(face = "bold")
  )

# Table summarizing the lactation dataset
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


model.colors <- c("black", "deepskyblue3", "darkseagreen4", 
                "coral3", "plum4", "darkorange3")

pars <- c("mu_alpha", "beta_pregnant_mod", "beta_lactating_mod")
bracket <- 
  list(
    c("Reproductive Effects on Viral Detection", 
      "Pregnancy Effect", 
      "Lactation Effect")
  )

# Generate a data frame of parameter means and HPDIs for each relevant
# parameter from each fit model

tidy.model.output <- data.frame(NULL)
conf.level <- 0.95

for (model.name in processed.model.names[1:6]) {
  
  temp.df <- 
    get(model.name) %>%
    select(all_of(pars)) %>%
    # Get means and HPDIs
    tidyr::pivot_longer(cols = everything(), names_to = "term", values_to = "value") %>%
    group_by(term) %>%
    summarize(
      estimate = mean(value),
      conf.low = HPDI(value, prob = conf.level)[1],
      conf.high = HPDI(value, prob = conf.level)[2]
    ) %>%
    ungroup() %>%
    # Clean up model names
    mutate(
      model = rep(model.name, length(pars)),
      model = str_replace(model, "model\\.f", ""),
      model = str_replace(model, "\\.p", ""),
      model = str_replace(model, "\\.", ""),
      model = ifelse(model == "", "All Viral Families", model)
    )
  
  tidy.model.output <- rbind(tidy.model.output, temp.df)
}

helper.data <- dat.f.trim %>%
  group_by(test_requested_viral_family) %>%
  count() %>%
  bind_rows(
    ., 
    data.frame(
      test_requested_viral_family = "All Viral Families", 
      n = nrow(dat.f.trim)
    )
  )

tidy.model.output <- tidy.model.output %>%
  left_join(
    ., helper.data, 
    by = c("model" = "test_requested_viral_family")
  ) %>%
  mutate(model = paste0(model, " (n = ", n, ")"))

models.excluding.all <- 
  unique(tidy.model.output$model)[!str_detect(unique(tidy.model.output$model), "All Viral Families")]
all.model <- 
  unique(tidy.model.output$model)[str_detect(unique(tidy.model.output$model), "All Viral Families")]

tidy.model.output <- tidy.model.output %>%
  mutate(
    model = factor(model, levels = c(all.model, models.excluding.all))
  )


plot1 <- {
  tidy.model.output %>%
    dotwhisker::relabel_predictors(
      c(mu_alpha = "Global Intercept",
        beta_pregnant_mod = "Pregnancy Effect",
        beta_lactating_mod = "Lactation Effect"
      )
    ) %>%
    dotwhisker::dwplot(
      dot_args = list(size = 2),
    ) +
    custom_theme +
    xlab("Parameter Estimate") +
    theme(
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9)
    ) +
    geom_vline(xintercept = 0, colour = "black", linetype = 2) +
    xlim(-10, 5) +
    scale_color_manual(
      values = rev(model.colors),
      name = "Viral Dataset"
    )
  } %>%
  dotwhisker::add_brackets(bracket)

plot2 <- model.f.p %>%
  select(matches("alpha_host_species\\[")) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "term", values_to = "value") %>%
  group_by(term) %>%
  summarize(
    estimate = mean(value),
    conf.low = HPDI(value, prob = conf.level)[1],
    conf.high = HPDI(value, prob = conf.level)[2]
  ) %>%
  ungroup() %>%
  dotwhisker::dwplot(
    .,
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

model.f.p %>%
  select(matches("alpha_host_species\\[")) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "term", values_to = "value") %>%
  group_by(term) %>%
  summarize(
    estimate = mean(value),
    conf.low = HPDI(value, prob = conf.level)[1],
    conf.high = HPDI(value, prob = conf.level)[2]
  ) %>%
  ungroup() %>%
  mutate(implied_probability = logistic(estimate)) %>%
  arrange(implied_probability)


plot_grid(plot1, plot2, 
          ncol = 1, scale = c(1, 0.95), rel_heights = c(1, 1.1),
          labels = "auto", label_size = 24)

ggsave("outputs/pooled_effects_model/Fig2.png", 
       width = 8, height = 10, dpi = 350)

plot1

ggsave("outputs/pooled_effects_model/Fig2.png", 
       width = 8, height = 6, dpi = 350)

# What proportion of posterior probability mass supports a negative pregnancy
# effect for each model?
sum(model.f.p$beta_pregnant_mod < 0)/
  nrow(model.f.p)
sum(model.f.Adenoviridae.p$beta_pregnant_mod < 0)/
  nrow(model.f.Adenoviridae.p)
sum(model.f.Coronaviridae.p$beta_pregnant_mod < 0)/
  nrow(model.f.Coronaviridae.p)
sum(model.f.Herpesviridae.p$beta_pregnant_mod < 0)/
  nrow(model.f.Herpesviridae.p)
sum(model.f.Paramyxoviridae.p$beta_pregnant_mod < 0)/
  nrow(model.f.Paramyxoviridae.p)
sum(model.f.Polyomaviridae.p$beta_pregnant_mod < 0)/
  nrow(model.f.Polyomaviridae.p)


# Replicate analysis and figures for varying slopes model

pars <- c("beta[1]", "beta[2]", "beta[3]")
bracket <- 
  list(
    c("Reproductive Effects on Viral Detection", 
      "Pregnancy Effect", 
      "Lactation Effect")
  )

tidy.model.output <- data.frame(NULL)
conf.level <- 0.95

for (model.name in processed.model.names[7:12]) {
  
  temp.df <-
    get(model.name) %>%
    select(all_of(pars)) %>%
    # Get means and HPDIs
    tidyr::pivot_longer(cols = everything(), names_to = "term", values_to = "value") %>%
    group_by(term) %>%
    summarize(
      estimate = mean(value),
      conf.low = HPDI(value, prob = conf.level)[1],
      conf.high = HPDI(value, prob = conf.level)[2]
    ) %>%
    ungroup() %>%
    # Clean up model names
    mutate(
      model = rep(model.name, length(pars)),
      model = str_replace(model, "model\\.f.v", ""),
      model = str_replace(model, "\\.p", ""),
      model = str_replace(model, "\\.", ""),
      model = ifelse(model == "", "All Viral Families", model)
    )
  
  tidy.model.output <- rbind(tidy.model.output, temp.df)
}

tidy.model.output <- tidy.model.output %>%
  left_join(
    ., helper.data, 
    by = c("model" = "test_requested_viral_family")
  ) %>%
  mutate(model = paste0(model, " (n = ", n, ")"))

models.excluding.all <- 
  unique(tidy.model.output$model)[!str_detect(unique(tidy.model.output$model), "All Viral Families")]
all.model <- 
  unique(tidy.model.output$model)[str_detect(unique(tidy.model.output$model), "All Viral Families")]

tidy.model.output <- tidy.model.output %>%
  mutate(
    model = factor(model, levels = c(all.model, models.excluding.all))
  )


plot1 <- {
  tidy.model.output %>%
    dotwhisker::relabel_predictors(
      c(`beta[1]` = "Intercept",
        `beta[2]` = "Pregnancy Effect",
        `beta[3]` = "Lactation Effect"
      )
    ) %>%
    dotwhisker::dwplot(
      dot_args = list(size = 2),
    ) +
    custom_theme +
    xlab("Parameter Estimate") +
    theme(
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9)
    ) +
    geom_vline(xintercept = 0, colour = "black", linetype = 2) +
    xlim(-10, 5) +
    scale_color_manual(
      values = rev(model.colors),
      name = "Viral Dataset"
    )
} %>%
  dotwhisker::add_brackets(bracket)

plot2 <- model.f.v.p %>%
  select(matches("beta_host_species\\[1")) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "term", values_to = "value") %>%
  group_by(term) %>%
  summarize(
    estimate = mean(value),
    conf.low = HPDI(value, prob = conf.level)[1],
    conf.high = HPDI(value, prob = conf.level)[2]
  ) %>%
  ungroup() %>%
  dotwhisker::dwplot(
    .,
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

model.f.v.p %>%
  select(matches("beta_host_species\\[1")) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "term", values_to = "value") %>%
  group_by(term) %>%
  summarize(
    estimate = mean(value),
    conf.low = HPDI(value, prob = conf.level)[1],
    conf.high = HPDI(value, prob = conf.level)[2]
  ) %>%
  ungroup() %>%
  mutate(implied_probability = logistic(estimate)) %>%
  arrange(implied_probability)


plot_grid(plot1, plot2, 
          ncol = 1, scale = c(1, 0.95), rel_heights = c(1, 1.1),
          labels = "auto", label_size = 24)

ggsave("outputs/varying_slopes_model/Fig2.png", 
       width = 8, height = 10, dpi = 350)

plot1

ggsave("outputs/varying_slopes_model/Fig2.png", 
       width = 8, height = 6, dpi = 350)

# What proportion of posterior probability mass supports a negative pregnancy
# effect for each model?
sum(model.f.v.p$`beta[2]` < 0)/
  nrow(model.f.v.p)
sum(model.f.v.Adenoviridae.p$`beta[2]` < 0)/
  nrow(model.f.v.Adenoviridae.p)
sum(model.f.v.Coronaviridae.p$`beta[2]` < 0)/
  nrow(model.f.v.Coronaviridae.p)
sum(model.f.v.Herpesviridae.p$`beta[2]` < 0)/
  nrow(model.f.v.Herpesviridae.p)
sum(model.f.v.Paramyxoviridae.p$`beta[2]` < 0)/
  nrow(model.f.v.Paramyxoviridae.p)
sum(model.f.v.Polyomaviridae.p$`beta[2]` < 0)/
  nrow(model.f.v.Polyomaviridae.p)

# What proportion of posterior probability mass supports a negative lactation
# effect for each model?
sum(model.f.v.p$`beta[3]` < 0)/
  nrow(model.f.v.p)
sum(model.f.v.Adenoviridae.p$`beta[3]` < 0)/
  nrow(model.f.v.Adenoviridae.p)
sum(model.f.v.Coronaviridae.p$`beta[3]` < 0)/
  nrow(model.f.v.Coronaviridae.p)
sum(model.f.v.Herpesviridae.p$`beta[3]` < 0)/
  nrow(model.f.v.Herpesviridae.p)
sum(model.f.v.Paramyxoviridae.p$`beta[3]` < 0)/
  nrow(model.f.v.Paramyxoviridae.p)
sum(model.f.v.Polyomaviridae.p$`beta[3]` < 0)/
  nrow(model.f.v.Polyomaviridae.p)
 
# /*
#==============================================================================
# */


# Figure 3 - Model-based predictions


sims <- 50
samples.per.sim <- 1000


big.sim.df <- data.frame(NULL)

for(x in processed.model.names[1:6]) {
  
  dat <- get(x)
  
  sim.df <- data.frame(
    rep(x, sims*3),
    rep(c("Non-reproductive", "Pregnant", "Lactating"), sims),
    rep(NA, sims*3),
    rep(NA, sims*3)
  )
  colnames(sim.df) <- c("model", "condition", "positives", "detection_prob")
  
  set.seed(8)
  
  for (i in (1:sims)*3 - 2) {
    
    sim.df$positives[i] <-
      sum(
        rbinom(samples.per.sim, 
               prob = sample(logistic(dat$mu_alpha), samples.per.sim), 
               size = 1)
      )
    
    sim.df$positives[i + 1] <-
      sum(
        rbinom(samples.per.sim, 
               prob = sample(logistic(dat$mu_alpha + dat$beta_pregnant_mod), 
                             samples.per.sim), 
               size = 1)
      )
    
    sim.df$positives[i + 2] <-
      sum(
        rbinom(samples.per.sim, 
               prob = sample(logistic(dat$mu_alpha + dat$beta_lactating_mod),
                             samples.per.sim), 
               size = 1)
      )
    
    sim.df$detection_prob[i:(i + 2)] <- sim.df$positives[i:(i + 2)]/samples.per.sim
  }
  
  big.sim.df <- bind_rows(big.sim.df, sim.df)
}

big.sim.df <- big.sim.df %>%
  mutate(
    model = 
      fct_relevel(model, "model.f.p"),
    condition = 
      fct_relevel(condition, c("Non-reproductive", "Pregnant", "Lactating"))
  )

model.labels <- c(
  model.f.p = "All Viral Families",
  model.f.Adenoviridae.p = "Adenoviridae",
  model.f.Coronaviridae.p = "Coronaviridae",
  model.f.Herpesviridae.p = "Herpesviridae",
  model.f.Paramyxoviridae.p = "Paramyxoviridae",
  model.f.Polyomaviridae.p = "Polyomaviridae"
)

alpha <- 1
color.values <- c(
  alpha("cornsilk3", alpha),
  alpha("darkseagreen", alpha),
  alpha("lightsteelblue", alpha)
)

plot <- big.sim.df %>%
  ggplot(aes(x = condition, y = detection_prob, color = condition)) +
  geom_jitter(height = 0) +
  ylab("Predicted Viral Detection Probability") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
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
  facet_wrap(~model, labeller = labeller(model = model.labels), nrow = 1)
  
summary.data <- big.sim.df %>%
  group_by(model, condition) %>%
  summarize(detection_prob = mean(detection_prob)) %>%
  mutate(x = rep(1:3)) %>%
  ungroup()


plot + 
  geom_path(
    data = summary.data, inherit.aes = FALSE,
    aes(x = x, y = detection_prob),
    color = alpha("black", 0.8), 
    linewidth = 0.7,
    arrow = arrow(angle = 20, length = unit(0.1, "inches"), type = "closed")
  )

ggsave("outputs/pooled_effects_model/Fig3.png", 
       height = 5, width = 8, dpi = 350)


# Replicate analysis and figure for varying slopes model

big.sim.df <- data.frame(NULL)

for(x in processed.model.names[7:12]) {
  
  dat <- get(x)
  
  sim.df <- data.frame(
    rep(x, sims*3),
    rep(c("Non-reproductive", "Pregnant", "Lactating"), sims),
    rep(NA, sims*3),
    rep(NA, sims*3)
  )
  colnames(sim.df) <- c("model", "condition", "positives", "detection_prob")
  
  set.seed(8)
  
  for (i in (1:sims)*3 - 2) {
    
    sim.df$positives[i] <-
      sum(
        rbinom(samples.per.sim, 
               prob = sample(logistic(dat$`beta[1]`), samples.per.sim), 
               size = 1)
      )
    
    sim.df$positives[i + 1] <-
      sum(
        rbinom(samples.per.sim, 
               prob = sample(logistic(dat$`beta[1]` + dat$`beta[2]`), 
                             samples.per.sim), 
               size = 1)
      )
    
    sim.df$positives[i + 2] <-
      sum(
        rbinom(samples.per.sim, 
               prob = sample(logistic(dat$`beta[1]` + dat$`beta[3]`),
                             samples.per.sim), 
               size = 1)
      )
    
    sim.df$detection_prob[i:(i + 2)] <- sim.df$positives[i:(i + 2)]/samples.per.sim
  }
  
  big.sim.df <- bind_rows(big.sim.df, sim.df)
}

big.sim.df <- big.sim.df %>%
  mutate(
    model = 
      fct_relevel(model, "model.f.v.p"),
    condition = 
      fct_relevel(condition, c("Non-reproductive", "Pregnant", "Lactating"))
  )

model.labels <- c(
  model.f.v.p = "All Viral Families",
  model.f.v.Adenoviridae.p = "Adenoviridae",
  model.f.v.Coronaviridae.p = "Coronaviridae",
  model.f.v.Herpesviridae.p = "Herpesviridae",
  model.f.v.Paramyxoviridae.p = "Paramyxoviridae",
  model.f.v.Polyomaviridae.p = "Polyomaviridae"
)

plot <- big.sim.df %>%
  ggplot(aes(x = condition, y = detection_prob, color = condition)) +
  geom_jitter(height = 0) +
  ylab("Predicted Viral Detection Probability") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
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
  facet_wrap(~model, labeller = labeller(model = model.labels), nrow = 1)

summary.data <- big.sim.df %>%
  group_by(model, condition) %>%
  summarize(detection_prob = mean(detection_prob)) %>%
  mutate(x = rep(1:3)) %>%
  ungroup()


plot + 
  geom_path(
    data = summary.data, inherit.aes = FALSE,
    aes(x = x, y = detection_prob),
    color = alpha("black", 0.8), 
    linewidth = 0.7,
    arrow = arrow(angle = 20, length = unit(0.1, "inches"), type = "closed")
  )

ggsave("outputs/varying_slopes_model/Fig3.png", 
       height = 5, width = 8, dpi = 350)


# Alternate plot for varying slopes model

big.df <- data.frame(NULL)

for(x in processed.model.names[7:12]) {
  
  dat <- get(x)
  
  temp.df <- data.frame(
    rep(x, nrow(dat)*3),
    rep(c("Non-reproductive", "Pregnant", "Lactating"), each = nrow(dat))
  )
  colnames(temp.df) <- c("model", "condition")
  
  temp.df$value <- c(
    logistic(dat$`beta[1]`),
    logistic(dat$`beta[1]` + dat$`beta[2]`), 
    logistic(dat$`beta[1]` + dat$`beta[3]`)
  )
  
  temp.df <- temp.df %>%
    group_by(model, condition) %>%
    summarize(
      mean = mean(value),
      median = median(value),
      mode = chainmode(value),
      lower_50 = HPDI(value, prob = 0.5)[1],
      upper_50 = HPDI(value, prob = 0.5)[2],
      lower_70 = HPDI(value, prob = 0.7)[1],
      upper_70 = HPDI(value, prob = 0.7)[2],
      lower_80 = HPDI(value, prob = 0.8)[1],
      upper_80 = HPDI(value, prob = 0.8)[2],
      lower_90 = HPDI(value, prob = 0.9)[1],
      upper_90 = HPDI(value, prob = 0.9)[2],
      lower_99 = HPDI(value, prob = 0.99)[1],
      upper_99 = HPDI(value, prob = 0.99)[2]
    ) %>%
    ungroup()

  big.df <- bind_rows(big.df, temp.df)
}

big.df <- big.df %>%
  mutate(
    model = 
      fct_relevel(model, "model.f.v.p"),
    condition = 
      fct_relevel(condition, c("Non-reproductive", "Pregnant", "Lactating"))
  ) %>%
  arrange(model, condition) %>%
  mutate(x = rep(1:3, times = 6))

big.df %>%
  ggplot(aes(x = condition, y = mode, color = condition)) +
  geom_point(size = 3) +
  geom_path(
    inherit.aes = FALSE, 
    aes(x = x, y = mode),
    color = "black", 
    linewidth = 0.3
  ) +
  geom_point(size = 3) +
  geom_linerange(aes(ymin = lower_80, ymax = upper_80), linewidth = 0.5) +
  geom_linerange(aes(ymin = lower_50, ymax = upper_50), linewidth = 1.5) +
  ylab("Estimated Viral Detection Probability") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
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
  facet_wrap(~model, labeller = labeller(model = model.labels), nrow = 1)

ggsave("outputs/varying_slopes_model/Fig3.png", 
       height = 5, width = 8, dpi = 350)

# /*
#==============================================================================
# */


# Supplementary Table 1


temp_table <- dat.f.trim %>%
  group_by(binomial) %>%
  summarize(
    tests = n(),
    positive_tests = sum(virus_detected == 1),
    tests_perc = positive_tests/tests,
    pregnant_tests = sum(pregnant_mod == 1),
    positive_pregnant_tests = sum(pregnant_mod == 1 & virus_detected == 1),
    pregnant_tests_perc = positive_pregnant_tests/pregnant_tests,
    lactating_tests = sum(lactating_mod == 1),
    positive_lactating_tests = sum(lactating_mod == 1 & virus_detected == 1),
    lactating_tests_perc = positive_lactating_tests/lactating_tests,
  )

totals_table <- data.frame(
  tests = sum(temp_table$tests),
  positive_tests = sum(temp_table$positive_tests),
  pregnant_tests = sum(temp_table$pregnant_tests),
  positive_pregnant_tests = sum(temp_table$positive_pregnant_tests),
  lactating_tests = sum(temp_table$lactating_tests),
  positive_lactating_tests = sum(temp_table$positive_lactating_tests)
) %>%
  mutate(
    tests_perc = positive_tests/tests,
    pregnant_tests_perc = positive_pregnant_tests/pregnant_tests,
    lactating_tests_perc = positive_lactating_tests/lactating_tests,
  ) %>%
  mutate_if(is.numeric, round, digits = 2) %>%
  mutate_all(as.character) %>%
  mutate(
    tests_perc = str_replace(tests_perc, "NaN", "NA"),
    pregnant_tests_perc = str_replace(pregnant_tests_perc, "NaN", "NA"),
    lactating_tests_perc = str_replace(lactating_tests_perc, "NaN", "NA"),
    tests = paste0(positive_tests, " / ", tests, " (", tests_perc, ")"),
    pregnant_tests = paste0(positive_pregnant_tests, " / ", pregnant_tests, " (", pregnant_tests_perc, ")"),
    lactating_tests = paste0(positive_lactating_tests, " / ", lactating_tests, " (", lactating_tests_perc, ")")
  ) %>%
  mutate(binomial = "Total") %>%
  select(binomial, tests, pregnant_tests, lactating_tests) %>%
  dplyr::rename(
    "Host Species" = binomial,
    "All Viral Tests" = tests,
    "Viral Tests from Pregnant Individuals" = pregnant_tests,
    "Viral Tests from Lactating Individuals" = lactating_tests)
  
temp_table <- temp_table %>%
  ungroup() %>%
  mutate_if(is.numeric, round, digits = 2) %>%
  mutate_all(as.character) %>%
  mutate(
    tests_perc = str_replace(tests_perc, "NaN", "NA"),
    pregnant_tests_perc = str_replace(pregnant_tests_perc, "NaN", "NA"),
    lactating_tests_perc = str_replace(lactating_tests_perc, "NaN", "NA"),
    tests = paste0(positive_tests, " / ", tests, " (", tests_perc, ")"),
    pregnant_tests = paste0(positive_pregnant_tests, " / ", pregnant_tests, " (", pregnant_tests_perc, ")"),
    lactating_tests = paste0(positive_lactating_tests, " / ", lactating_tests, " (", lactating_tests_perc, ")")
  ) %>%
  select(binomial, tests, pregnant_tests, lactating_tests) %>%
  dplyr::rename(
    "Host Species" = binomial,
    "All Viral Tests" = tests,
    "Viral Tests from Pregnant Individuals" = pregnant_tests,
    "Viral Tests from Lactating Individuals" = lactating_tests) %>%
  arrange(`Host Species`)

temp_table %>%
  bind_rows(totals_table) %>%
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

par(
  mar = c(0, 0, 0, 0) + 0.05,
  family = "Arial"
)

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
  tibble(
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
  tibble(
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


p <- bayesplot::mcmc_trace(
  model.f %>% 
    posterior::as_draws(),
  pars = c(
    "mu_alpha", "beta_pregnant_mod", "beta_lactating_mod",
    "sigma_vector[1]", "sigma_vector[2]", "sigma_vector[3]",
    "sigma_vector[4]", "sigma_vector[5]", "sigma_vector[6]"
  )
)

levels(p$data$parameter) <- c(
  "Global Intercept", 
  "Pregnancy Effect", "Lactation Effect",
  "σ (Host Species)", "σ (Year)", "σ (Country)", 
  "σ (Specimen Type)", "σ (Viral Test Protocol)", "σ (Diagnostic Laboratory)"
)


p + 
  scale_color_viridis(discrete = TRUE, option = "plasma", end = 0.85) +
  ggtitle("Parameter trace plots for Bayesian model of viral detection in adult female bats") +
  theme(
    text = element_text(size = 14, color = "black", family = "sans"),
    plot.title = element_text(size = 18),
    strip.text.x = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

ggsave("outputs/pooled_effects_model/FigS4.png", 
       width = 10, height = 10, dpi = 350)


# Replicate for the varying slopes model

p <- bayesplot::mcmc_trace(
  model.f.v %>% 
    posterior::as_draws(),
  pars = c(
    "beta[1]", "beta[2]", "beta[3]",
    "sigma_host_species[1]", "sigma_host_species[2]", "sigma_host_species[3]",
    "sigma_vector[1]", "sigma_vector[2]", "sigma_vector[3]",
    "sigma_vector[4]", "sigma_vector[5]"
  ),
  facet_args = list(ncol = 3)
)

levels(p$data$parameter) <- c(
  "Community Intercept", 
  "Community Pregnancy Effect", "Community Lactation Effect",
  "σ (Host Species Intercepts)", "σ (Host Species Pregnancy Effects)", 
  "σ (Host Species Lactation Effects)",
  "σ (Year)", "σ (Country)", "σ (Specimen Type)", "σ (Viral Test Protocol)", 
  "σ (Diagnostic Laboratory)"
)


p + 
  scale_color_viridis(discrete = TRUE, option = "plasma", end = 0.85) +
  ggtitle("Parameter trace plots for Bayesian model of viral detection in adult female bats") +
  theme(
    text = element_text(size = 14, color = "black", family = "sans"),
    plot.title = element_text(size = 18),
    strip.text.x = element_text(size = 14, face = "bold"),
    legend.position = "none"
  )

ggsave("outputs/varying_slopes_model/FigS4.png", 
       width = 12, height = 12, dpi = 350)

# /*
#==============================================================================
# */


# Figure S6 - Ridgeline plots for varying effects


dat.df <- dat.f.trim
model.df <- model.f.p

cols.to.plot <- grep("alpha", colnames(model.df), value = T) %>%
  grep("mu_|tilde|species\\[", ., value = T, invert = T)

var.effect.labels <- 
  list(
    "host_species_offset" = sort(unique(dat.df$binomial)),
    "year" = levels(factor(dat.df$year)),
    "country" = levels(factor(dat.df$country)),
    "specimen_type_group" = sort(unique(dat.df$specimen_type_group)),
    "test_requested_mod" = sort(unique(dat.df$test_requested_mod)),
    "diagnostic_laboratory_name" = sort(unique(dat.df$diagnostic_laboratory_name))
  )

label.df <- read_csv("data/lookup_tables/laboratory_name_cleanup.csv")

indexes <- 
  match(var.effect.labels$diagnostic_laboratory_name, label.df$diagnostic_laboratory_name)

var.effect.labels$diagnostic_laboratory_name <- 
  label.df$diagnostic_laboratory_name_mod[indexes]

assert_that(length(cols.to.plot) == length(flatten(var.effect.labels)))


# To plot all varying effects

png("outputs/pooled_effects_model/FigS6a.png", width = 1200, height = 800)

p <- model.df %>% 
  select(all_of(cols.to.plot)) %>%
  gather("parameter", "value", cols.to.plot, factor_key = T) %>%
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
  xlim(-11, 11) +
  # To facet by varying effects group
  facet_wrap(~varying_effect_group, scale = "free_y", 
             labeller = as_labeller(var.effect.group.labels)) +
  theme_ridges(center_axis_labels = TRUE, font_size = 16) +
  # To remove y-axis labels
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_blank(), 
    axis.title = element_text(size = 32),
    plot.tag = element_text(size = 50, face = "bold"),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 18, face = "bold", margin = margin(b = 5))
  ) +
  coord_cartesian(clip = "off") +
  labs(tag = "a")

p

dev.off()


# To plot varying effects one at a time

plotting.list <- list(
  c("host_species_offset", "year", "country",
    "specimen_type_group", "test_requested_mod", "diagnostic_laboratory_name"),
  c("outputs/pooled_effects_model/FigS6b.png", "outputs/pooled_effects_model/FigS6c.png",
    "outputs/pooled_effects_model/FigS6d.png", "outputs/pooled_effects_model/FigS6e.png", 
    "outputs/pooled_effects_model/FigS6f.png", "outputs/pooled_effects_model/FigS6g.png"),
  c("Host Species", "Year of Sample Collection",
    "Country of Sample Collection", "Specimen Type", 
    "Viral Test Protocol", "Diagnostic Laboratory Conducting Testing"
     ),
  c(-6, -6, -6, -6, -11, -6)
)
labs <- c("b", "c", "d", "e", "f", "g")

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
      select(all_of(cols.to.plot)) %>%
      gather("parameter", "value", cols.to.plot, factor_key = T) %>%
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
      theme(
        axis.title = element_text(size = 32),
        plot.tag = element_text(size = 50, face = "bold"),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 16)
      ) +
      scale_y_discrete(expand = expansion(add = c(0.2, 1.5))) +
      labs(tag = labs[i])
  )
  
  dev.off()
}


# Replicate for the varying slopes model

dat.df <- dat.f.trim
model.df <- model.f.v.p

cols.to.plot <- grep("beta_host_species", colnames(model.df), value = T)

labels <- rep(sort(unique(dat.df$binomial)), each = 3)
names(labels) <- cols.to.plot

varying.intercept.slope.group.labels <- c(
  "beta_host_species[1" = "Species-Specific Intercepts",
  "beta_host_species[2" = "Species-Specific Pregnancy Effects",
  "beta_host_species[3" = "Species-Specific Lactation Effects"
) 


# To plot all varying intercepts and slopes by species

png("outputs/varying_slopes_model/FigS6.png", width = 1500, height = 1500)

p <- model.df %>% 
  select(all_of(cols.to.plot)) %>%
  gather("parameter", "value", cols.to.plot, factor_key = T) %>%
  mutate(
    varying_effect_group = gsub(",[0-9]+\\]", "", parameter),
    varying_effect_group = factor(varying_effect_group, names(varying.intercept.slope.group.labels))
  ) %>%
  # To order the distributions by median values
  mutate(parameter = reorder(parameter, value, median)) %>%
  mutate(host_species = plyr::mapvalues(parameter, cols.to.plot, labels)) %>%
  
  ggplot(aes(x = value, y = parameter, 
             height = ..density.., fill = varying_effect_group)) +
  xlab("Parameter Value") + ylab("") +
  geom_density_ridges() +
  scale_y_discrete(labels = labels) +
  # To facet by varying effects group
  facet_wrap(
    ~varying_effect_group, scale = "free", 
    labeller = labeller(
      varying_effect_group = varying.intercept.slope.group.labels
    )
  ) +
  theme_ridges(center_axis_labels = TRUE, font_size = 16) +
  # To remove y-axis labels
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(face = "italic"),
    axis.title = element_text(size = 32),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 18, face = "bold", margin = margin(b = 5))
  ) +
  coord_cartesian(clip = "off")

p

dev.off()


# To plot varying intercepts and slopes by species one at a time

plotting.list <- list(
  c("beta_host_species\\[1", "beta_host_species\\[2", "beta_host_species\\[3"),
  c("outputs/varying_slopes_model/FigS6a.png", "outputs/varying_slopes_model/FigS6b.png", 
    "outputs/varying_slopes_model/FigS6c.png"),
  c("Species-Specific Intercepts", "Species-Specific Pregnancy Effects",
    "Species-Specific Lactation Effects"),
  c(c(-12, 2), c(-4, 2), c(-5, 3))
)
labs <- c("a", "b", "c")

for (i in 1:length(plotting.list[[1]])) {
  
  cols.to.plot.sub <- grep(plotting.list[[1]][i], colnames(model.df), value = T)
  
  fill.color <- ggplot_build(p)$data[[1]] %>%
    distinct(fill) %>%
    slice(i) %>%
    unlist()
  
  png(plotting.list[[2]][i], width = 1000, height = 1500)
  
  print(
    model.df %>% 
      select(all_of(cols.to.plot)) %>%
      gather("parameter", "value", cols.to.plot, factor_key = T) %>%
      mutate(
        varying_effect_group = gsub(",[0-9]+\\]", "", parameter),
        varying_effect_group = factor(varying_effect_group, names(varying.intercept.slope.group.labels))
      ) %>%
      # To order the distributions by median values
      mutate(parameter = reorder(parameter, value, median)) %>%
      mutate(host_species = plyr::mapvalues(parameter, cols.to.plot, labels)) %>%
      filter(parameter %in% cols.to.plot.sub) %>%
      
      ggplot(aes(x = value, y = parameter, 
                 height = ..density..)) +
      xlab("Parameter Value") + ylab(plotting.list[[3]][i]) +
      geom_density_ridges(fill = fill.color) +
      coord_cartesian(xlim = c(plotting.list[[4]][i*2-1], plotting.list[[4]][i*2])) +
      theme_ridges(center_axis_labels = TRUE, font_size = 12) +
      theme(
        axis.title = element_text(size = 32),
        plot.tag = element_text(size = 50, face = "bold"),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 16, face = "italic")
      ) +
      scale_y_discrete(labels = labels, expand = expansion(add = c(0.2, 1.5))) +
      labs(tag = labs[i])
  )
  
  dev.off()
}


cols.to.plot <- grep("alpha", colnames(model.df), value = T) %>%
  grep("mu_|tilde|species\\[", ., value = T, invert = T)

var.effect.labels <- 
  list(
    "year" = levels(factor(dat.df$year)),
    "country" = levels(factor(dat.df$country)),
    "specimen_type_group" = sort(unique(dat.df$specimen_type_group)),
    "test_requested_mod" = sort(unique(dat.df$test_requested_mod)),
    "diagnostic_laboratory_name" = sort(unique(dat.df$diagnostic_laboratory_name))
  )

label.df <- read_csv("data/lookup_tables/laboratory_name_cleanup.csv")

indexes <- 
  match(var.effect.labels$diagnostic_laboratory_name, label.df$diagnostic_laboratory_name)

var.effect.labels$diagnostic_laboratory_name <- 
  label.df$diagnostic_laboratory_name_mod[indexes]

assert_that(length(cols.to.plot) == length(flatten(var.effect.labels)))


# To plot all other varying effects

png("outputs/varying_slopes_model/FigS7a.png", width = 1200, height = 800)

p <- model.df %>% 
  select(all_of(cols.to.plot)) %>%
  gather("parameter", "value", cols.to.plot, factor_key = T) %>%
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
  xlim(-11, 11) +
  # To facet by varying effects group
  facet_wrap(~varying_effect_group, scale = "free_y", 
             labeller = as_labeller(var.effect.group.labels)) +
  theme_ridges(center_axis_labels = TRUE, font_size = 16) +
  # To remove y-axis labels
  theme(
    axis.text.x = element_text(size = 20),
    axis.text.y = element_blank(), 
    axis.title = element_text(size = 32),
    plot.tag = element_text(size = 50, face = "bold"),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_text(size = 18, face = "bold", margin = margin(b = 5))
  ) +
  coord_cartesian(clip = "off") +
  labs(tag = "a")

p

dev.off()


# To plot other varying effects one at a time

plotting.list <- list(
  c("year", "country",
    "specimen_type_group", "test_requested_mod", "diagnostic_laboratory_name"),
  c("outputs/varying_slopes_model/FigS7b.png", "outputs/varying_slopes_model/FigS7c.png",
    "outputs/varying_slopes_model/FigS7d.png", "outputs/varying_slopes_model/FigS7e.png", 
    "outputs/varying_slopes_model/FigS7f.png"),
  c("Year of Sample Collection",
    "Country of Sample Collection", "Specimen Type", 
    "Viral Test Protocol", "Diagnostic Laboratory Conducting Testing"
  ),
  c(-6, -6, -6, -11, -6)
)
labs <- c("b", "c", "d", "e", "f")

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
      select(all_of(cols.to.plot)) %>%
      gather("parameter", "value", cols.to.plot, factor_key = T) %>%
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
      theme(
        axis.title = element_text(size = 32),
        plot.tag = element_text(size = 50, face = "bold"),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 16)
      ) +
      scale_y_discrete(expand = expansion(add = c(0.2, 1.5))) +
      labs(tag = labs[i])
  )
  
  dev.off()
}

# /*
#==============================================================================
# */


# Figure S5 - In-sample prediction plot

# Generate data frame of all alpha values from the full fit model
d.preds <- model.f %>%
  posterior::as_draws_df() %>%
  select(contains("alpha")) %>%
  select(!contains("_")) 

# How many iterations total?
n.iter <- dim(d.preds)[1]

# How many data points total?
n.datapoints <- dim(d.preds)[2]

# Pivot the alpha matrix such that we get a data frame with columns for
# the model iteration, data point, and probability of success
d.preds <- d.preds %>%
  tidyr::pivot_longer(
    cols = contains("alpha"),
    names_to = "data point",
    values_to = "alpha"
  ) %>%
  mutate(
    iteration = rep(1:n.iter, each = n.datapoints),
    prob = logistic(alpha)
  )

set.seed(1)

# Generate 0/1 predictions for each data point/iteration combination using
# the modeled probability of success
d.preds$pred <- rbinom(
  n = length(d.preds$prob), 
  size = 1, 
  prob = d.preds$prob
)

# Record the test requested viral family for each observation
d.preds$viral_family <- rep(
  dat.f.trim$test_requested_viral_family,
  times = n.iter
)

# Summarize the predictions (across all viral families) to generate metrics
# of test positivity
full.dataset.preds <- d.preds %>%
  group_by(iteration) %>%
  summarize(
    n = n(),
    n_positive = sum(pred),
    positivity = n_positive/n
  ) %>%
  ungroup() %>%
  mutate(viral_family = rep("All Data", times = n()))

# Summarize the predictions by viral family to generate metrics
# of test positivity and bind in the all viral families data
viral.family.preds <- d.preds %>%
  group_by(viral_family, iteration) %>%
  summarize(
    n = n(),
    n_positive = sum(pred),
    positivity = n_positive/n
  ) %>%
  ungroup() %>%
  bind_rows(full.dataset.preds)

# Summarize test positivity predictions using HPDIs
plot.data <- viral.family.preds %>%
  group_by(viral_family, n) %>%
  mutate(
    mean = mean(positivity),
    lower.95 = HPDI(positivity, prob = 0.95)[1],
    upper.95 = HPDI(positivity, prob = 0.95)[2],
    lower.50 = HPDI(positivity, prob = 0.5)[1],
    upper.50 = HPDI(positivity, prob = 0.5)[2]
  ) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    viral_family = paste0(viral_family, "\n(n = ", n, ")"),
    viral_family = fct_relevel(
      viral_family,
      "All Data\n(n = 9694)"
    )
  ) %>%
  select(-c(iteration, n, n_positive, positivity))

# Get analogous test positivity metrics for the observed data
full.dataset.obs <- dat.f.trim %>%
  group_by() %>%
  summarize(
    n = n(),
    n_positive = sum(virus_detected),
    positivity = n_positive/n
  ) %>%
  ungroup() %>%
  mutate(viral_family = rep("All Data", times = n()))

viral.family.obs <- dat.f.trim %>%
  group_by(test_requested_viral_family) %>%
  summarize(
    n = n(),
    n_positive = sum(virus_detected),
    positivity = n_positive/n
  ) %>%
  ungroup() %>%
  rename(viral_family = test_requested_viral_family) %>%
  bind_rows(full.dataset.obs) %>%
  mutate(
    viral_family = paste0(viral_family, "\n(n = ", n, ")"),
    viral_family = fct_relevel(
      viral_family,
      "All Data\n(n = 9694)"
    )
  )

plot.data <- plot.data %>%
  left_join(
    ., viral.family.obs,
    by = "viral_family"
  ) %>%
  mutate(
    in_50_interval = ifelse(
      positivity <= upper.50 & positivity >= lower.50,
      TRUE, FALSE
    )
  )

# Plot and save
plot.data %>%
  ggplot(aes(x = viral_family)) +
  geom_linerange(
    aes(ymin = lower.95, ymax = upper.95), 
    linewidth = 1, color = "darkgrey"
  ) +
  geom_linerange(
    aes(ymin = lower.50, ymax = upper.50), 
    linewidth = 4, color = "darkgrey"
  ) +
  geom_point(
    aes(y = positivity), 
    size = 2, color = "darkred"
  ) +
  geom_vline(xintercept = 1.5, lty = 2) +
  ylab("Test positivity") +
  xlab("") +
  ylim(0, 0.2) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(face = "bold")
  )

ggsave("outputs/pooled_effects_model/FigS5.png", 
       height = 5, width = 10, dpi = 350)


# Replicate for the varying slopes model

# Generate data frame of all alpha values from the full fit model
d.preds <- model.f.v %>%
  posterior::as_draws_df() %>%
  select(contains("alpha")) %>%
  select(!contains("_")) 

# How many iterations total?
n.iter <- dim(d.preds)[1]

# How many data points total?
n.datapoints <- dim(d.preds)[2]

# Pivot the alpha matrix such that we get a data frame with columns for
# the model iteration, data point, and probability of success
d.preds <- d.preds %>%
  tidyr::pivot_longer(
    cols = contains("alpha"),
    names_to = "data point",
    values_to = "alpha"
  ) %>%
  mutate(
    iteration = rep(1:n.iter, each = n.datapoints),
    prob = logistic(alpha)
  )

set.seed(1)

# Generate 0/1 predictions for each data point/iteration combination using
# the modeled probability of success
d.preds$pred <- rbinom(
  n = length(d.preds$prob), 
  size = 1, 
  prob = d.preds$prob
)

# Record the test requested viral family for each observation
d.preds$viral_family <- rep(
  dat.f.trim$test_requested_viral_family,
  times = n.iter
)

# Summarize the predictions (across all viral families) to generate metrics
# of test positivity
full.dataset.preds <- d.preds %>%
  group_by(iteration) %>%
  summarize(
    n = n(),
    n_positive = sum(pred),
    positivity = n_positive/n
  ) %>%
  ungroup() %>%
  mutate(viral_family = rep("All Data", times = n()))

# Summarize the predictions by viral family to generate metrics
# of test positivity and bind in the all viral families data
viral.family.preds <- d.preds %>%
  group_by(viral_family, iteration) %>%
  summarize(
    n = n(),
    n_positive = sum(pred),
    positivity = n_positive/n
  ) %>%
  ungroup() %>%
  bind_rows(full.dataset.preds)

# Summarize test positivity predictions using HPDIs
plot.data <- viral.family.preds %>%
  group_by(viral_family, n) %>%
  mutate(
    mean = mean(positivity),
    lower.95 = HPDI(positivity, prob = 0.95)[1],
    upper.95 = HPDI(positivity, prob = 0.95)[2],
    lower.50 = HPDI(positivity, prob = 0.5)[1],
    upper.50 = HPDI(positivity, prob = 0.5)[2]
  ) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    viral_family = paste0(viral_family, "\n(n = ", n, ")"),
    viral_family = fct_relevel(
      viral_family,
      "All Data\n(n = 9694)"
    )
  ) %>%
  select(-c(iteration, n, n_positive, positivity))

plot.data <- plot.data %>%
  left_join(
    ., viral.family.obs,
    by = "viral_family"
  ) %>%
  mutate(
    in_50_interval = ifelse(
      positivity <= upper.50 & positivity >= lower.50,
      TRUE, FALSE
    )
  )

# Plot and save
plot.data %>%
  ggplot(aes(x = viral_family)) +
  geom_linerange(
    aes(ymin = lower.95, ymax = upper.95), 
    linewidth = 1, color = "darkgrey"
  ) +
  geom_linerange(
    aes(ymin = lower.50, ymax = upper.50), 
    linewidth = 4, color = "darkgrey"
  ) +
  geom_point(
    aes(y = positivity),
    size = 2, color = "darkred"
  ) +
  geom_vline(xintercept = 1.5, lty = 2) +
  ylab("Test positivity") +
  xlab("") +
  ylim(0, 0.2) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(face = "bold")
  )

ggsave("outputs/varying_slopes_model/FigS5.png", 
       height = 5, width = 10, dpi = 350)
