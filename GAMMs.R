# This script contains a GAMM analysis of a visual world paradigm study
# on pronoun and anaphor comprehension, run with German-speaking children and adults.
# Created by Nevena Klobučar.
# Updated: 26.03.2025.

# Load packages
library(tidyverse)
library(mgcv)
library(itsadug)

# Set wd
setwd("C:/Users/neven/Dropbox/DPBE_German3_JB2_deploy/myfiles")

# Import eye tracking data
eyedata <- read.table("analysis.txt", 
                      sep = ",", 
                      header = TRUE, 
                      stringsAsFactors = TRUE, 
                      strip.white = TRUE, 
                      fill = TRUE, 
                      quote = "")


# ..............................................................................
#  Define dependent variable: Target vs. Competitor ----
# ..............................................................................

# Change TRUE/FALSE values in Aerea of Interest columns to 1/0
eyedata <- eyedata %>%
  mutate(
    Target = as.numeric(Target),     
    Competitor = as.numeric(Competitor) 
  )

# Convert data set from wider to longer
eyedata <- eyedata %>%
  pivot_longer(cols = c(Target, Competitor),  
               names_to = "IA",             
               values_to = "Fixation")   %>%
  arrange(IA, Time)

# Interest Area as factor and check levels
eyedata$IA <- factor(eyedata$IA, levels = c("Target", "Competitor"))
levels(eyedata$IA)

# Create a data set where only fixations to target are counted. 
# 1 for target means 0 for competitor.
targetdata <- eyedata %>% 
  
  # Keep only rows where either the target or competitor were fixated.
  filter(Fixation == 1) %>%
  
  # Create a new variable indicating whether the target was fixated.
  mutate(Target = ifelse(IA == "Target", 1, 0)) %>%
  
  # Drop unnecessary information.
  dplyr::select(-c(IA, Fixation)) %>%
  distinct(.keep_all = TRUE) %>%
  droplevels()

# ..............................................................................
#  We will analyze children's and adults data seperately ----
# ..............................................................................

# Create subset of adult data
adults_eyedata <- targetdata %>%
  filter(group == "adults")

#create subset of child data
children_eyedata <- targetdata %>%
  filter(group == "children")

# ..............................................................................
# First, we analyze adults data
# ...................................
# Logit-transform looks to Target 
# ..............................................................................

# Keep bins of 20 ms (50 Hz)
binned_data <- adults_eyedata %>%
  group_by(Participant, Items, Time, condition) %>% 
  summarise(
    looks = sum(Target),
    total = n(),
    .groups = "drop"
  )

# Compute logit
c <- 0.5  # small constant to avoid log(0)
binned_data <- binned_data %>%
  mutate(
    elogit = log((looks + c) / ((total - looks) + c)),
    weight = 1 / (looks + c) + 1 / ((total - looks) + c)
  )

# ..............................................................................
# Convert categorical predictors to factors ----
# ..............................................................................

binned_data <- binned_data %>%
  mutate(
    condition = as.factor(condition),
    Participant = as.factor(Participant),
    Items = as.factor(Items)
  )

# order data (necessary for checking autocorrelation in residuals)
binned_data <- as.data.frame(binned_data)

binned_data <- droplevels(binned_data[order(binned_data$Participant, 
                                            binned_data$Items, 
                                            binned_data$Time),])

summary(binned_data)

# ..............................................................................
# Define window of interest ----
# ..............................................................................

# The window of interest includes time from -500 to 1100 ms (end of sentence) from pronoun or anaphor onset
binned_data <- droplevels(binned_data[binned_data$Time > -500 & binned_data$Time < 1100,])

# ..............................................................................
# Check for autocorrelation and define rho ----
# ..............................................................................

# Create start event

binned_data <- binned_data %>%
  mutate(EventID = interaction(Participant, Items, drop = TRUE))
str(binned_data$EventID)


binned_data <- start_event(
  binned_data,
  column = "Time",
  event = "EventID",
  label.event = "Event"
)

# Run plain model
m0 <- bam(elogit ~ s(Time) + #smooth for time across conditions
            s(Participant, bs = 're') + #random effects for participants
            s(Items, bs = 're'), #random effects for items
          data = binned_data)

# Check for autocorrelation
acf(resid(m0), main="acf(resid(m0))")


# Determine rho parameter
r1 <- start_value_rho(m0, plot=TRUE)

# r1=0.94

# ..............................................................................
# Run GAMMs ----
# ..............................................................................

# Check levels of condition (anaphor_true as reference)
binned_data$condition <- relevel(binned_data$condition, ref="anaphor_true")
# inspect contrasts:
contrasts(binned_data$condition)

m1 <- bam(
  elogit ~ condition + 
    s(Time, by = condition)  + #separate smooths over time for each condition
    s(Participant, bs = 're') + #random smooths for participants
    s(Items, bs = 're'), #random smooths for items
  data = binned_data,
  rho = r1,
  AR.start = start.event # handle temporal autocorrelation in residuals
)

summary(m1)

# plot (transform elogit back to probability)
invlogit <- function(x) 1 / (1 + exp(-x))

plot_smooth(
  m1,
  view = "Time",
  plot_all = "condition",
  transform = invlogit,
  rug = FALSE,
  main = "Proportion of looks to target",
  ylab = "Proportion of target looks"
)

# pairwise comparison of true conditions
plot_diff(m1, 
          view = "Time", 
          comp = list(condition = c("anaphor_true", "pronoun_true")), 
          main = "Anaphor vs Pronoun (True)")

# pairwise comparison of true conditions
plot_diff(
  m1,
  view = "Time",
  comp = list(condition = c("anaphor_false", "pronoun_false")),
  main = "Anaphor vs Pronoun (False)"
)


# ..............................................................................
#  We now analyze children's data  ----
# ..............................................................................

# Repeat all previous steps

# Keep bins of 20 ms (50 Hz)
binned_children <- children_eyedata %>%
  group_by(Participant, Items, Time, condition) %>% 
  summarise(
    looks = sum(Target),
    total = n(),
    .groups = "drop"
  )

# Compute logit
c <- 0.5  # small constant to avoid log(0)
binned_children <- binned_children %>%
  mutate(
    elogit = log((looks + c) / ((total - looks) + c)),
    weight = 1 / (looks + c) + 1 / ((total - looks) + c)
  )

# ..............................................................................
# Convert categorical predictors to factors ----
# ..............................................................................

binned_children <- binned_children %>%
  mutate(
    condition = as.factor(condition),
    Participant = as.factor(Participant),
    Items = as.factor(Items)
  )

# ..............................................................................
# Define window of interest ----
# ..............................................................................

# The window of interest includes time from -500 to 1100 ms (end of sentence) from pronoun or anaphor onset
binned_children <- droplevels(binned_children[binned_children$Time > -500 & binned_children$Time < 1100,])


# order data (necessary for checking autocorrelation in residuals)
binned_children <- as.data.frame(binned_children)

binned_children <- droplevels(binned_children[order(binned_children$Participant, 
                                            binned_children$Items, 
                                            binned_children$Time),])

head(binned_children)


# ..............................................................................
# Check for autocorrelation and define rho ----
# ..............................................................................

# Create start event

binned_children <- binned_children %>%
  mutate(EventID = interaction(Participant, Items, drop = TRUE))
str(binned_children$EventID)


binned_children <- start_event(
  binned_children,
  column = "Time",
  event = "EventID",
  label.event = "Event"
)


# Run plain model
m0_children <- bam(elogit ~ s(Time) + #smooth for time across conditions
            s(Participant, bs = 're') + #random effects for participants
            s(Items, bs = 're'), #random effects for items
          data = binned_children)

# Check for autocorrelation
acf(resid(m0_children), main="acf(resid(m0_children))")


# Determine rho parameter
r1b <- start_value_rho(m0_children, plot=TRUE)

# r1b=0.94

# ..............................................................................
# Run GAMMs ----
# ..............................................................................

# Check levels of condition (anaphor_true as reference)
binned_children$condition <- relevel(binned_children$condition, ref="anaphor_true")
# inspect contrasts:
contrasts(binned_children$condition)

m1_children <- bam(
  elogit ~ condition + 
    s(Time, by = condition)  + #separate smooths over time for each condition
    s(Participant, bs = 're') + #random smooths for participants
    s(Items, bs = 're'), #random smooths for items
  data = binned_children,
  rho = r1,
  AR.start = start.event # handle temporal autocorrelation in residuals
)

summary(m1_children)

# plot (transform elogit back to probability)
invlogit <- function(x) 1 / (1 + exp(-x))

plot_smooth(
  m1_children,
  view = "Time",
  plot_all = "condition",
  transform = invlogit,
  rug = FALSE,
  main = "Proportion of looks to target",
  ylab = "Proportion of target looks"
)

# pairwise comparison of true conditions
plot_diff(m1_children, 
          view = "Time", 
          comp = list(condition = c("anaphor_true", "pronoun_true")), 
          main = "Anaphor vs Pronoun (True)")

# pairwise comparison of false conditions
plot_diff(
  m1_children,
  view = "Time",
  comp = list(condition = c("anaphor_false", "pronoun_false")),
  main = "Anaphor vs Pronoun (False)"
)
