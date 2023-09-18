#' ---
#' title: "Bat Viral Detection and Reproduction Analysis - Data Prep"
#' author: "EcoHealth Alliance M&A Team - Code Drafted by Evan Eskew"
#' ---

# /*
#==============================================================================
# */


#+ data_import_chunk, echo=FALSE, message=FALSE, results="hide"


# Load packages and functions

library(tidyverse)
library(assertthat)
library(reskew) # devtools::install_github("eveskew/reskew")

source("R/functions.R")


# Import raw P1 data tables

e <- read_csv("data/p1_extracts/events.csv") %>%
  remove_NA_cols()
a <- read_csv("data/p1_extracts/animals.csv") %>%
  remove_NA_cols()
s <- read_csv("data/p1_extracts/specimens.csv") %>%
  remove_NA_cols()
t <- read_csv("data/p1_extracts/tests.csv") %>%
  remove_NA_cols()
v <- read_csv("data/p1_extracts/viruses.csv") %>%
  remove_NA_cols()
ts <- read_csv("data/p1_extracts/test_specimen_ids.csv") %>%
  remove_NA_cols()


# Add higher-level bat taxonomic information to the animal table

a <- left_join(
  a, read_csv("data/lookup_tables/P1_bat_classification.csv"), 
  by = "family"
)

# Add information on viral family taxonomy to the test table

t <- left_join(
  t, read_csv("data/lookup_tables/P1_viral_family_of_tests.csv"),
  by = "test_requested"
)

# Simplify viral species names in a new column in the viruses table
# Note: this is because anything called "strain of..." should
# not be considered a distinct viral species

v <- v %>% 
  mutate(viral_species = stringi::stri_replace_first_regex(
    virus_name, "(new\\s)?strain\\sof\\s", ""))


# Create full P1 data frame by joining individual tables

d <- full_join(e, a, by = "event_id") %>%
  full_join(s, by = "animal_id") %>%
  full_join(ts, by = "specimen_id") %>%
  full_join(t, by = "test_id") %>%
  full_join(v, by = "test_id")

# Purely aesthetic label changes

d <- d %>%
  mutate(
    specimen_type_group = case_when(
      specimen_type_group == "oral/nasal/oropharyngeal swabs" ~ 
        "oral/nasal/oropharyngeal swab",
      TRUE ~ specimen_type_group
    )
  )


# Create data frame to keep track of sample sizes for intermediate datasets

d.sample.sizes <- data.frame(
  "d", 
  n_distinct(d$animal_id, na.rm = T), 
  n_distinct(d$binomial, na.rm = T)
)

colnames(d.sample.sizes) <- c("data_frame", "n_animals", "n_species")

d.sample.sizes$data_frame <- as.character(d.sample.sizes$data_frame)

# /*
#==============================================================================
# */


#+ data_modification_chunk, echo=FALSE, results="hide"


d2 <- d


# Replace erroneous "sample_date" entries for relevant sites

d2$sample_date[which(d2$site_name == "BR_DF_3DEC2012_FAZENDAUFAM_site1_rainy" &
                       d2$sample_date == "2012-03-12")] <- "2012-12-03"

d2$sample_date[which(d2$site_name == "BR_DF_9DEC2012_FAZENDAUFAM_site3_rainy" &
                       d2$sample_date == "2012-09-12")] <- "2012-12-09"

# Replace erroneous "specimen_date" entries of "1900-01-01" with NAs

d2$specimen_date[which(d2$specimen_date == "1900-01-01")] <- NA


# Create a modified sampling date column that is equal to "sample_date"
# when it exists and equal to "event_date" in other cases

d2$date_mod <-
  ifelse(is.na(d2$sample_date),
         as.character(d2$event_date), as.character(d2$sample_date))

# Substitute "specimen_date" in "date_mod" if it exists and 
# "sample_date" doesn't

d2$date_mod <- 
  ifelse((is.na(d2$sample_date) & !is.na(d2$specimen_date)),
         as.character(d2$specimen_date), as.character(d2$date_mod))


# Add more granular date info

d2$day_of_year <- yday(d2$date_mod)
d2$month <- month(d2$date_mod) %>%
  as.factor()
d2$year <- year(d2$date_mod)


# Modify country names

d2 <- d2 %>%
  mutate(
    country = case_when(
      country == "Malaysia, Peninsular" ~ "Malaysia",
      country == "Malaysia, Sabah" ~ "Malaysia",
      TRUE ~ country
    )
  )


# Modify reproductive trait variables to resolve ambiguities

d2 <- d2 %>%
  mutate(
    pregnant_mod = case_when(
      pregnant == "N/A" ~ NA_real_,
      pregnant == "Unknown" ~ NA_real_,
      pregnant == "Yes" ~ 1,
      pregnant == "No" ~ 0
    ),
    lactating_mod = case_when(
      lactating == "N/A" ~ NA_real_,
      lactating == "Unknown" ~ NA_real_,
      lactating == "Yes" ~ 1,
      lactating == "No" ~ 0
    )
  )

table(d2$pregnant_mod, useNA = "ifany")
table(d2$lactating_mod, useNA = "ifany")


# Modify "test_requested" column to incorporate distinct testing protocols

d2 <- d2 %>%
  mutate(
    test_requested_mod = 
      paste0(test_requested, ": ",  test_requested_protocol),
  ) %>%
  left_join(
    ., read_csv("data/lookup_tables/test_requested_mod_cleanup.csv"),
    by = "test_requested_mod"
  ) %>%
  mutate(
    test_requested_mod = 
      ifelse(
        !is.na(test_requested_mod_new), 
        test_requested_mod_new, 
        test_requested_mod
      )
  )

# Create the "confirmation_result_mod" column, using the "confirmation_result" 
# to determine positives

d2$confirmation_result_mod <- 
  ifelse(d2$confirmation_result == "Positive", "Positive", "Negative")

# Convert "confirmation_result_mod" into a binary numeric variable

d2$virus_detected <- ifelse(d2$confirmation_result_mod == "Positive", 1, 0)


# Modify misspellings of viral family names

d2 <- d2 %>%
  mutate(
    viral_family = case_when(
      viral_family == "Flaviridae" ~ "Flaviviridae",
      viral_family == "Togavirdae" ~ "Togaviridae",
      TRUE ~ viral_family
    )
  )

table(d2$viral_family)

assert_that(
  sum(
    filter(d2, !is.na(confirmation_result) & !is.na(viral_family)) %>%
      pull(test_requested_viral_family) != 
      filter(d2, !is.na(confirmation_result) & !is.na(viral_family)) %>%
      pull(viral_family)
  ) == 0
)


# Add to sample size data frame

d.sample.sizes <- 
  rbind(d.sample.sizes, 
        c("d2", n_distinct(d2$animal_id, na.rm = T), 
          n_distinct(d2$binomial, na.rm = T)))

# /*
#==============================================================================
# */


#+ data_filtering_chunk, echo=FALSE, results="hide"


# Filter data to bats

d.bat <- d2 %>%
  filter(order == "Chiroptera") %>%
  droplevels()

nrow(d.bat)


# Add to sample size data frame

d.sample.sizes <- 
  rbind(d.sample.sizes, 
        c("d.bat", n_distinct(d.bat$animal_id, na.rm = T), 
          n_distinct(d.bat$binomial, na.rm = T)))


# Create further filtered data frame for analysis

# Make sure PREDICT protocol == "TRUE"
# Make sure animal_classification == "Wild"
# Remove age_class == NA
# Remove age_class == "Unknown"
# Remove serology, next generation sequencing, and real time PCR tests
# Remove test_requested_protocol == "other"
# Remove pooled tests
# Remove confirmation_result == NA
# Remove confirmation_result == "Pool positive - do not count"

d3 <- d.bat %>%
  filter(
    predict_protocol == 1,
    animal_classification == "Wild",
    !is.na(age_class),
    age_class != "Unknown",
    test_type_broad != "Serologic_Tests",
    test_type_broad != "Sequencing",
    test_type_specific != "Real time PCR",
    test_requested_protocol != "other",
    pooled != 1,
    !is.na(confirmation_result),
    confirmation_result != "Pool positive - do not count"
  ) %>%
  droplevels()

nrow(d3)
table(d3$predict_protocol, useNA = "ifany")
table(d3$animal_classification, useNA = "ifany")
table(d3$age_class, useNA = "ifany")
table(d3$test_type_broad, useNA = "ifany")
table(d3$test_type_specific, useNA = "ifany")
table(d3$test_result, useNA = "ifany")
table(d3$confirmation_result, useNA = "ifany")
table(d3$confirmation_result_mod, useNA = "ifany")
table(d3$pooled, useNA = "ifany")


# Since the same "test_id" can produce multiple "viral_species" hits,
# only take the first "viral_species" for a given "test_id" so that each
# row represents a single test

d3 <- d3 %>%
  arrange(desc(test_id), desc(viral_species)) %>%
  distinct(test_id, .keep_all = TRUE)

nrow(d3)
assert_that(nrow(d3) == nrow(distinct(d3, test_id)))


# Add to sample size data frame

d.sample.sizes <- 
  rbind(d.sample.sizes, 
        c("d3", n_distinct(d3$animal_id, na.rm = T), 
          n_distinct(d3$binomial, na.rm = T)))

# /*
#==============================================================================
# */


#+ data_definition_chunk, echo=FALSE, results="hide"


# Adult female data

dat.f <- d3 %>%
  # Filter to adult females
  filter(age_class %in% c("Subadult", "Adult") & sex == "Female") %>%
  # Make sure all records have reproductive and species identity data
  filter(
    !is.na(pregnant_mod),
    !is.na(lactating_mod),
    !is.na(binomial)
  ) %>%
  arrange(animal_id, specimen_id, test_requested, test_requested_protocol)

# Find tests that have no viral detections for filtering

tests.no.positives <- dat.f %>%
  group_by(test_requested) %>%
  summarize(positives = sum(virus_detected)) %>%
  filter(positives == 0) %>%
  select(test_requested) %>%
  unlist(use.names = FALSE)

dat.f <- dat.f %>%
  filter(!(test_requested %in% tests.no.positives)) %>%
  droplevels()

nrow(dat.f)

table(dat.f$sex, dat.f$age_class, useNA = "ifany")
table(dat.f$age_class, dat.f$virus_detected, useNA = "ifany")
table(dat.f$pregnant_mod, useNA = "ifany")
table(dat.f$lactating_mod, useNA = "ifany")


# Save cleaned data

saveRDS(dat.f, file = "data/cleaned_data/dat.f.rds")

trim.cols <- c(
  "year", "country", 
  "animal_id", "order", "family", "genus", "binomial", 
  "pregnant_mod", "lactating_mod", 
  "specimen_id", "specimen_type_group",
  "test_requested", "test_requested_protocol", 
  "test_requested_mod", "test_requested_viral_family",
  "diagnostic_laboratory_name", "test_id", 
  "confirmation_result_mod", "virus_detected", "viral_species"
)

dat.f.trim <- select(dat.f, all_of(trim.cols))

saveRDS(dat.f.trim, file = "data/cleaned_data/dat.f.trim.rds")
write_csv(dat.f.trim, file = "data/cleaned_data/dat.f.trim.csv")


# Add to sample size data frame

d.sample.sizes <- 
  rbind(d.sample.sizes, 
        c("dat.f.trim", n_distinct(dat.f.trim$animal_id, na.rm = T), 
          n_distinct(dat.f.trim$binomial, na.rm = T)))


# Save sample size data frame

write_csv(d.sample.sizes, file = "outputs/d.sample.sizes.csv")


# Create viral family-specific data subsets

viral.fam.preg.pos <- dat.f.trim %>% 
  group_by(test_requested_viral_family, pregnant_mod) %>%
  summarize(positive = sum(virus_detected)) %>%
  filter(pregnant_mod == 1, positive >= 1) %>%
  pull(test_requested_viral_family)

viral.fam.lac.pos <- dat.f.trim %>% 
  group_by(test_requested_viral_family, lactating_mod) %>%
  summarize(positive = sum(virus_detected)) %>%
  filter(lactating_mod == 1, positive >= 1) %>%
  pull(test_requested_viral_family)

viral.families <- intersect(viral.fam.preg.pos, viral.fam.lac.pos)

data.list <- vector("list", length(viral.families) + 1)
data.list[[1]] <- dat.f.trim
names(data.list)[1] <- "dat.f"

for(i in seq_along(viral.families)) {
  
  assign(paste0("dat.f.", viral.families[i]), 
         filter(dat.f.trim, test_requested_viral_family == viral.families[i]) %>%
           droplevels()
  )
  
  data.list[[i + 1]] <- get(paste0("dat.f.", viral.families[i]))
  
  names(data.list)[i + 1] <- paste0("dat.f.", viral.families[i])
}

# /*
#==============================================================================
# */


#+ data_definition_chunk_stan, echo=FALSE, results="hide"


# Define Stan data for all data frames

for(i in seq_along(data.list)) {
  
  saveRDS(get_stan_data(data.list[[i]]), 
          file = paste0("stan/cleaned_data/", names(data.list)[i], ".stan.rds")
  )
}

# /*
#==============================================================================
# */


#+ data_summary_chunk, echo=FALSE, results="hide"


# Summarize adult female dataset

summarize(dat.f.trim, 
          sample_size = n(), 
          n_viral_families = 
            n_distinct(test_requested_viral_family, na.rm = TRUE),
          n_viruses = n_distinct(viral_species, na.rm = TRUE),
          n_species = n_distinct(binomial, na.rm = TRUE),
          n_years = n_distinct(year, na.rm = TRUE),
          n_countries = n_distinct(country, na.rm = TRUE),
          n_specimen_types = n_distinct(specimen_type_group, na.rm = TRUE),
          n_test_protocols = n_distinct(test_requested_mod, na.rm = TRUE),
          n_diagnostic_labs = n_distinct(diagnostic_laboratory_name, 
                                         na.rm = TRUE)
)
