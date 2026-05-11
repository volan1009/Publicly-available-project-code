# ==============================================================================
# Methodological Pilot: Single-Outcome GBTM Stress Test
# Author: Wondering Data Ghost #2
# Objective: Simulate a recurrent psychiatric count outcome under observation
#            loss and administrative censoring, then evaluate class recovery.
# ==============================================================================
# install.packages("flexmix")

library(dplyr)
library(tidyr)
library(flexmix)
library(ggplot2)

# 1. Define pilot simulation parameters ----------------------------------------
set.seed(202605)             # Fix the seed so the pilot is reproducible.
N <- 5000                    # Use a smaller cohort than the planned full study.
T_max <- 10                  # Follow each subject for up to 10 intervals.
classes <- 3                 # Fit three latent trajectory classes.
pi_k <- c(0.5, 0.3, 0.2)     # Set the true class proportions.

# Use class-specific quadratic Poisson trajectories for the true outcome process.
beta_matrix <- data.frame(
  class = 1:3,
  beta0 = c(-1.5, -0.5, 0.5),
  beta1 = c(0.1, 0.3, 0.2),
  beta2 = c(-0.01, -0.05, -0.02)
)

# 2. Generate the fully observed ground truth ----------------------------------
generate_ground_truth <- function(N, T_max, pi_k, beta_matrix) {
  # Sample one latent class per subject.
  latent_class <- sample(1:length(pi_k), size = N, replace = TRUE, prob = pi_k)

  # Build the person-period grid and attach the class-specific coefficients.
  df_grid <- expand_grid(id = 1:N, t = 1:T_max) %>%
    left_join(data.frame(id = 1:N, class = latent_class), by = "id") %>%
    left_join(beta_matrix, by = "class")

  # Draw the true recurrent count outcome from the class-specific Poisson mean.
  df_truth <- df_grid %>%
    mutate(
      log_lambda = beta0 + beta1 * t + beta2 * (t^2),
      lambda = exp(log_lambda),
      y_true = rpois(n(), lambda)
    ) %>%
    select(id, true_class = class, t, y_true)

  return(df_truth)
}

df_truth <- generate_ground_truth(N, T_max, pi_k, beta_matrix)

# 3. Apply observation loss and administrative censoring -----------------------
apply_data_corruption <- function(df_truth, T_max, trigger_threshold = 3) {
  # Approximate broad observation loss with one simple cutoff rule.
  df_after_loss <- df_truth %>%
    group_by(id) %>%
    mutate(
      cum_observed_events = cumsum(y_true),
      observation_loss_flag = lag(cum_observed_events, default = 0) >= trigger_threshold,
      y_after_observation_loss = ifelse(observation_loss_flag, 0L, y_true)
    ) %>%
    ungroup()

  # Use a simple right-skewed pilot distribution for follow-up truncation.
  subject_ids <- sort(unique(df_truth$id))
  censor_times <- pmin(T_max + 2, 4 + rnbinom(length(subject_ids), size = 2, mu = 2))

  # Treat all values after the censoring time as missing rather than zero.
  df_corrupted <- df_after_loss %>%
    left_join(data.frame(id = subject_ids, C_i = censor_times), by = "id") %>%
    mutate(
      y_observed = ifelse(t > C_i, NA, y_after_observation_loss)
    ) %>%
    select(id, true_class, t, y_true, y_observed)

  return(df_corrupted)
}

df_corrupted <- apply_data_corruption(df_truth, T_max)

# Inspect the first few rows of the final observed pilot dataset.
print(head(df_corrupted, 10))

# 4. Define helper functions for fitting and evaluation -------------------------
fit_poisson_gbtm <- function(df, outcome_name, classes) {
  # Fit the same quadratic Poisson mixture form described in the manuscript.
  flexmix(
    as.formula(paste(outcome_name, "~ t + I(t^2) | id")),
    data = df,
    k = classes,
    model = FLXMRglm(family = "poisson"),
    control = list(iter.max = 500)
  )
}

extract_id_assignments <- function(df, model) {
  # Collapse the row-level labels back to one class assignment per subject.
  df %>%
    mutate(raw_pred_class = clusters(model)) %>%
    distinct(id, raw_pred_class)
}

make_relabel_map <- function(true_class, pred_class, n_classes) {
  # Apply a greedy one-to-one relabeling rule based on maximum overlap.
  overlap_table <- table(
    True = factor(true_class, levels = 1:n_classes),
    Pred = factor(pred_class, levels = 1:n_classes)
  )

  relabel_map <- c()
  temp_table <- overlap_table

  for (i in seq_len(n_classes)) {
    max_idx <- which(temp_table == max(temp_table), arr.ind = TRUE)[1, ]
    true_c <- as.numeric(rownames(temp_table)[max_idx[1]])
    pred_c <- as.numeric(colnames(temp_table)[max_idx[2]])

    relabel_map[as.character(pred_c)] <- true_c
    temp_table[max_idx[1], ] <- -1
    temp_table[, max_idx[2]] <- -1
  }

  return(relabel_map)
}

extract_trajectories <- function(model, relabel_map, type_label, T_max) {
  # Recover the fitted Poisson mean curve from each mixture component.
  coefs <- parameters(model)
  rows <- vector("list", ncol(coefs))

  # Locate the intercept, linear, and quadratic terms in the parameter matrix.
  idx_b0 <- grep("Intercept", rownames(coefs))
  idx_b1 <- grep("t$", rownames(coefs))
  idx_b2 <- grep("t\\^2", rownames(coefs))

  for (comp_idx in seq_len(ncol(coefs))) {
    # Map the raw component label back to the aligned class label.
    true_c <- as.numeric(relabel_map[as.character(comp_idx)])
    t_seq <- 1:T_max

    # Convert the fitted linear predictor into the expected Poisson count.
    lambda <- exp(
      coefs[idx_b0, comp_idx] +
        coefs[idx_b1, comp_idx] * t_seq +
        coefs[idx_b2, comp_idx] * (t_seq^2)
    )

    rows[[comp_idx]] <- data.frame(
      t = t_seq,
      true_class = true_c,
      Trajectory_Type = type_label,
      Expected_Count = lambda
    )
  }

  bind_rows(rows)
}

# 5. Fit the pilot GBTM to the degraded observed data --------------------------
print("--- Fitting the pilot GBTM to degraded observed data ---")

# Drop administratively censored values before model fitting.
df_fit <- df_corrupted %>%
  filter(!is.na(y_observed))

# Fit the three-class quadratic Poisson model to the observed counts.
set.seed(202605)
gbtm_model <- fit_poisson_gbtm(df_fit, "y_observed", classes)

# Convert the grouped fit into one recovered class per subject.
df_pred <- extract_id_assignments(df_fit, gbtm_model)

# Join the recovered classes back to the known ground truth.
df_eval <- df_corrupted %>%
  distinct(id, true_class) %>%
  left_join(df_pred, by = "id") %>%
  drop_na(raw_pred_class)

# Relabel the recovered classes so they align with the true class numbering.
relabel_map <- make_relabel_map(df_eval$true_class, df_eval$raw_pred_class, classes)

df_eval <- df_eval %>%
  mutate(aligned_pred_class = as.numeric(relabel_map[as.character(raw_pred_class)]))

# Summarize the class-specific attrition rate used in the manuscript.
attrition_results <- df_eval %>%
  group_by(true_class) %>%
  summarise(
    total_in_true_class = n(),
    correctly_classified = sum(aligned_pred_class == true_class),
    attrition_rate = round(1 - (correctly_classified / total_in_true_class), 4),
    .groups = "drop"
  ) %>%
  mutate(
    class_label = case_when(
      true_class == 1 ~ "Low Burden",
      true_class == 2 ~ "Intermediate Burden",
      true_class == 3 ~ "High Burden"
    )
  )

print("--- Class-specific attrition rate after observation loss and censoring ---")
print(attrition_results)

# Inspect how one true class is redistributed after the data are degraded.
print(
  df_eval %>%
    filter(true_class == 2) %>%
    count(aligned_pred_class) %>%
    mutate(prop = n / sum(n))
)

# 6. Fit the same model to complete data for comparison ------------------------
print("--- Fitting the same pilot GBTM to complete data ---")

# Fit the model to the fully observed truth to create a benchmark comparison.
set.seed(202605)
perfect_model <- fit_poisson_gbtm(df_truth, "y_true", classes)

# Extract one recovered class label per subject from the complete-data fit.
df_pred_complete <- extract_id_assignments(df_truth, perfect_model)

# Build the evaluation frame for the complete-data relabeling step.
df_eval_complete <- df_truth %>%
  distinct(id, true_class) %>%
  left_join(df_pred_complete, by = "id")

# Relabel the complete-data fit using the same overlap-based rule.
relabel_map_complete <- make_relabel_map(
  df_eval_complete$true_class,
  df_eval_complete$raw_pred_class,
  classes
)

# 7. Visualize theoretical and fitted trajectories -----------------------------
# Compute the true expected count under the data-generating equation.
df_theoretical <- expand_grid(t = 1:T_max, true_class = 1:classes) %>%
  left_join(beta_matrix, by = c("true_class" = "class")) %>%
  mutate(
    Expected_Count = exp(beta0 + beta1 * t + beta2 * (t^2)),
    Trajectory_Type = "1. Theoretical truth"
  ) %>%
  select(t, true_class, Trajectory_Type, Expected_Count)

# Extract the fitted mean curves from the complete-data and degraded-data models.
df_complete_fit <- extract_trajectories(
  perfect_model,
  relabel_map_complete,
  "2. Fitted on complete data",
  T_max
)

df_degraded_fit <- extract_trajectories(
  gbtm_model,
  relabel_map,
  "3. Fitted after observation loss and censoring",
  T_max
)

# Stack the three trajectory sources into one plotting frame.
df_plot <- bind_rows(df_theoretical, df_complete_fit, df_degraded_fit) %>%
  mutate(
    class_label = factor(
      true_class,
      levels = 1:3,
      labels = c("Low Burden", "Intermediate Burden", "High Burden")
    )
  )

# Plot the theoretical curves against the two fitted recovery curves.
p <- ggplot(
  df_plot,
  aes(
    x = t,
    y = Expected_Count,
    color = Trajectory_Type,
    linetype = Trajectory_Type
  )
) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ class_label, scales = "fixed") +
  scale_color_manual(values = c("black", "#0072B2", "#D55E00")) +
  scale_linetype_manual(values = c("dashed", "solid", "solid")) +
  scale_x_continuous(breaks = 1:10) +
  labs(
    title = "Trajectory recovery under observation loss and administrative censoring",
    subtitle = "Theoretical truth versus fits on complete and degraded observed data",
    x = "Follow-up interval (t)",
    y = "Expected count",
    color = "Series",
    linetype = "Series"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "grey90"),
    panel.grid.minor = element_blank()
  )

print(p)
