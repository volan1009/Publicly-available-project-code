# Step 1. Check required packages.
required_pkgs <- c("gemtc", "coda", "ggplot2", "tidyr")

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop("Install required packages first: ", paste(missing_pkgs, collapse = ", "))
}

library(gemtc)
library(coda)
library(ggplot2)
library(tidyr)


# Step 2. Define path and output helpers.
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)

  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)))
  }

  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = TRUE)))
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

normalize_names <- function(x) {
  tolower(gsub("[^[:alnum:]]+", "", x))
}

slugify <- function(x) {
  gsub("(^_+|_+$)", "", tolower(gsub("[^[:alnum:]]+", "_", x)))
}

save_png <- function(file, plot_fun, width = 1800, height = 1200, res = 180) {
  png(file, width = width, height = height, res = res)
  on.exit(dev.off(), add = TRUE)
  plot_fun()
}

write_text_output <- function(expr, file) {
  capture.output(expr, file = file)
}

has_multi_arm_trials <- function(data) {
  any(table(data$study) > 2)
}

save_forest_plot <- function(relative_effects, file, digits = 3) {
  n_rows <- nvar(as.mcmc.list(relative_effects)[[1]])
  plot_height <- max(1600, 110 * n_rows + 500)

  save_png(
    file,
    function() forest(relative_effects, digits = digits, ask = FALSE),
    height = plot_height
  )
}

find_windows_jags_homes <- function() {
  roots <- unique(c(
    Sys.getenv("JAGS_HOME", unset = ""),
    file.path(Sys.getenv("ProgramFiles", unset = ""), "JAGS"),
    file.path(Sys.getenv("ProgramFiles(x86)", unset = ""), "JAGS")
  ))

  roots <- roots[nzchar(roots)]
  candidates <- character()

  for (root in roots) {
    if (!dir.exists(root)) {
      next
    }

    if (basename(root) == "JAGS") {
      subdirs <- list.dirs(root, recursive = FALSE, full.names = TRUE)
      candidates <- c(candidates, subdirs[grepl("^JAGS-4", basename(subdirs), ignore.case = TRUE)])
    } else {
      candidates <- c(candidates, root)
    }
  }

  unique(normalizePath(candidates[dir.exists(candidates)], winslash = "/", mustWork = FALSE))
}

prepend_to_path <- function(paths) {
  existing_paths <- paths[dir.exists(paths)]
  if (length(existing_paths) == 0) {
    return(invisible(NULL))
  }

  existing_paths <- normalizePath(existing_paths, winslash = "/", mustWork = FALSE)
  current_path <- strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]
  current_path <- current_path[nzchar(current_path)]

  Sys.setenv(PATH = paste(unique(c(existing_paths, current_path)), collapse = .Platform$path.sep))
}

ensure_rjags_ready <- function() {
  installed_packages <- rownames(installed.packages())

  if (!("rjags" %in% installed_packages)) {
    stop(
      "The R package 'rjags' is not installed.\n",
      "Install JAGS 4.x first, then run install.packages('rjags')."
    )
  }

  load_result <- tryCatch(loadNamespace("rjags"), error = function(e) e)
  if (!inherits(load_result, "error")) {
    return(invisible(TRUE))
  }

  if (.Platform$OS.type == "windows") {
    jags_homes <- find_windows_jags_homes()

    for (jags_home in jags_homes) {
      Sys.setenv(JAGS_HOME = jags_home)

      prepend_to_path(c(
        jags_home,
        file.path(jags_home, "bin"),
        file.path(jags_home, "x64", "bin"),
        file.path(jags_home, "i386", "bin")
      ))

      load_result <- tryCatch(loadNamespace("rjags"), error = function(e) e)
      if (!inherits(load_result, "error")) {
        return(invisible(TRUE))
      }
    }
  }

  stop(
    "rjags is installed but JAGS 4.x is not available.\n",
    "Install JAGS 4.x from https://sourceforge.net/projects/mcmc-jags/files/ and rerun.\n",
    "If JAGS is already installed, set JAGS_HOME to the JAGS installation folder before starting R.\n",
    "Original rjags error: ", conditionMessage(load_result)
  )
}


# Step 3. Read and validate arm-level input data.
validate_common_fields <- function(data) {
  data$study <- trimws(data$study)
  data$treatment <- tolower(trimws(data$treatment))

  if (any(data$study == "")) {
    stop("Blank study values detected.")
  }

  if (any(data$treatment == "")) {
    stop("Blank treatment values detected.")
  }

  if (any(duplicated(data[c("study", "treatment")]))) {
    stop("Duplicate study-treatment rows detected.")
  }

  invalid_ids <- unique(data$treatment[!grepl("^[a-z0-9_]+$", data$treatment)])
  if (length(invalid_ids) > 0) {
    stop("Invalid treatment IDs: ", paste(invalid_ids, collapse = ", "))
  }

  if (!"placebo" %in% data$treatment) {
    stop("Reference treatment 'placebo' not found.")
  }

  arm_counts <- table(data$study)
  if (any(arm_counts < 2)) {
    stop(
      "One-arm studies detected: ",
      paste(names(arm_counts[arm_counts < 2]), collapse = ", ")
    )
  }

  data
}

read_continuous_outcome <- function(path) {
  raw <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(raw) <- normalize_names(names(raw))

  required <- c("study", "treatment", "mean", "stddev", "samplesize")
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0) {
    stop("Missing required columns in ", basename(path), ": ", paste(missing, collapse = ", "))
  }

  data <- validate_common_fields(raw[required])
  data$mean <- as.numeric(data$mean)
  data$std.dev <- as.numeric(data$stddev)
  data$sampleSize <- as.integer(data$samplesize)
  data$stddev <- NULL
  data$samplesize <- NULL

  if (anyNA(data$mean) || anyNA(data$std.dev) || anyNA(data$sampleSize)) {
    stop("Continuous outcome contains NA after numeric coercion.")
  }

  data[c("study", "treatment", "mean", "std.dev", "sampleSize")]
}

read_binary_outcome <- function(path) {
  raw <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(raw) <- normalize_names(names(raw))

  required <- c("study", "treatment", "responders", "samplesize")
  missing <- setdiff(required, names(raw))
  if (length(missing) > 0) {
    stop("Missing required columns in ", basename(path), ": ", paste(missing, collapse = ", "))
  }

  data <- validate_common_fields(raw[required])
  data$responders <- as.integer(data$responders)
  data$sampleSize <- as.integer(data$samplesize)
  data$samplesize <- NULL

  if (anyNA(data$responders) || anyNA(data$sampleSize)) {
    stop("Binary outcome contains NA after numeric coercion.")
  }

  if (any(data$responders < 0) || any(data$responders > data$sampleSize)) {
    stop("Responder counts must be between 0 and sample size.")
  }

  data[c("study", "treatment", "responders", "sampleSize")]
}


# Step 4. Save cumulative ranking plots.
save_cumulative_rank_plot <- function(ranks, outcome_label, file) {
  ranks_matrix <- unclass(ranks)
  ranks_df <- data.frame(
    treatment = rownames(ranks_matrix),
    ranks_matrix,
    check.names = FALSE,
    row.names = NULL
  )

  ranks_long <- pivot_longer(
    ranks_df,
    cols = -treatment,
    names_to = "rank",
    values_to = "probability"
  )

  ranks_long$rank <- as.integer(gsub("[^0-9]+", "", ranks_long$rank))
  ranks_long <- ranks_long[order(ranks_long$treatment, ranks_long$rank), ]

  ranks_long$cumprob <- ave(
    ranks_long$probability,
    ranks_long$treatment,
    FUN = cumsum
  )

  p <- ggplot(ranks_long, aes(x = rank, y = cumprob, colour = treatment, group = treatment)) +
    geom_line(linewidth = 0.9) +
    labs(
      x = "Rank",
      y = "Cumulative probability",
      title = paste("Cumulative ranking plot:", outcome_label)
    ) +
    theme_bw()

  ggsave(filename = file, plot = p, width = 8, height = 5, dpi = 300)
}


# Step 5. Run the continuous-outcome workflow.
Sys.setenv(JAGS_HOME = "D:/Program Files/JAGS/JAGS-4.3.2")

run_continuous_nma <- function(data_file, outcome_label, preferred_direction = -1,
                               reference_treatment = "placebo", n_adapt = 20000,
                               n_iter = 50000, thin = 1, n_chain = 4, seed = 20240511) {
  ensure_rjags_ready()

  slug <- slugify(outcome_label)
  data <- read_continuous_outcome(data_file)
  network <- mtc.network(data.ab = data, description = outcome_label)

  # Save the basic network figure.
  save_png(file.path("figures", paste0("network_", slug, ".png")), function() plot(network))

  # Fit the consistency model and the UME model.
  model_consistency <- mtc.model(
    network,
    type = "consistency",
    n.chain = n_chain,
    likelihood = "normal",
    link = "identity",
    linearModel = "random"
  )

  set.seed(seed)
  results_consistency <- mtc.run(
    model_consistency,
    n.adapt = n_adapt,
    n.iter = n_iter,
    thin = thin
  )

  if (has_multi_arm_trials(data)) {
    results_ume <- NULL
    message("Skipping UME model for ", outcome_label, ": gemtc warns that UME does not handle multi-arm trials correctly.")
  } else {
    model_ume <- mtc.model(
      network,
      type = "ume",
      n.chain = n_chain,
      likelihood = "normal",
      link = "identity",
      linearModel = "random"
    )

    set.seed(seed)
    results_ume <- mtc.run(
      model_ume,
      n.adapt = n_adapt,
      n.iter = n_iter,
      thin = thin
    )
  }

  # Export model summaries and core diagnostics.
  write_text_output(summary(results_consistency), file.path("results", paste0("summary_", slug, ".txt")))
  if (!is.null(results_ume)) {
    write_text_output(summary(results_ume), file.path("results", paste0("summary_ume_", slug, ".txt")))
  } else {
    writeLines(
      c(
        "UME model not run.",
        "Reason: gemtc warns that UME does not handle multi-arm trials correctly for this network."
      ),
      file.path("results", paste0("summary_ume_", slug, ".txt"))
    )
  }

  relative_effects <- relative.effect(
    results_consistency,
    reference_treatment,
    preserve.extra = FALSE
  )

  save_forest_plot(
    relative_effects,
    file.path("figures", paste0("forest_", slug, ".png")),
    digits = 3
  )

  save_png(
    file.path("figures", paste0("trace_", slug, ".png")),
    function() plot(results_consistency)
  )

  save_png(
    file.path("figures", paste0("gelman_", slug, ".png")),
    function() gelman.plot(results_consistency)
  )

  # Export ranking results.
  ranks <- rank.probability(results_consistency, preferredDirection = preferred_direction)

  save_png(
    file.path("figures", paste0("rank_probability_", slug, ".png")),
    function() plot(ranks, beside = TRUE)
  )

  save_cumulative_rank_plot(ranks, outcome_label, file.path("figures", paste0("sucra_", slug, ".png")))

  sucra_values <- sucra(ranks)
  sucra_results <- data.frame(
    treatment = names(sucra_values),
    sucra = as.numeric(sucra_values),
    row.names = NULL
  )

  write.csv(sucra_results, file.path("results", paste0("sucra_", slug, ".csv")), row.names = FALSE)

  league_table <- round(relative.effect.table(results_consistency), 2)
  write.csv(league_table, file.path("results", paste0("league_table_", slug, ".csv")))

  # Run node-splitting when the network structure allows it.
  node_split <- tryCatch(
    mtc.nodesplit(
      network,
      n.adapt = n_adapt,
      n.iter = n_iter,
      thin = thin,
      n.chain = n_chain,
      likelihood = "normal",
      link = "identity",
      linearModel = "random"
    ),
    error = function(e) {
      message("Node-splitting skipped for ", outcome_label, ": ", e$message)
      NULL
    }
  )

  if (!is.null(node_split)) {
    node_summary <- summary(node_split)
    write_text_output(node_summary, file.path("results", paste0("nodesplit_", slug, ".txt")))

    save_png(
      file.path("figures", paste0("nodesplit_", slug, ".png")),
      function() plot(node_summary, digits = 5)
    )
  }

  # Summarize heterogeneity.
  heterogeneity <- mtc.anohe(
    network,
    n.adapt = n_adapt,
    n.iter = n_iter,
    thin = thin,
    n.chain = n_chain,
    likelihood = "normal",
    link = "identity",
    linearModel = "random"
  )

  heterogeneity_summary <- summary(heterogeneity)
  write_text_output(heterogeneity_summary, file.path("results", paste0("heterogeneity_", slug, ".txt")))

  save_png(
    file.path("figures", paste0("heterogeneity_", slug, ".png")),
    function() plot(heterogeneity_summary, digits = 5)
  )
}


# Step 6. Run the binary-outcome workflow.
run_binary_nma <- function(data_file, outcome_label, preferred_direction = -1,
                           reference_treatment = "placebo", n_adapt = 20000,
                           n_iter = 50000, thin = 1, n_chain = 4, seed = 20240511) {
  ensure_rjags_ready()

  slug <- slugify(outcome_label)
  data <- read_binary_outcome(data_file)
  network <- mtc.network(data.ab = data, description = outcome_label)

  # Save the basic network figure.
  save_png(file.path("figures", paste0("network_", slug, ".png")), function() plot(network))

  # Fit the consistency model and the UME model.
  model_consistency <- mtc.model(
    network,
    type = "consistency",
    n.chain = n_chain,
    likelihood = "binom",
    link = "logit",
    linearModel = "random"
  )

  set.seed(seed)
  results_consistency <- mtc.run(
    model_consistency,
    n.adapt = n_adapt,
    n.iter = n_iter,
    thin = thin
  )

  if (has_multi_arm_trials(data)) {
    results_ume <- NULL
    message("Skipping UME model for ", outcome_label, ": gemtc warns that UME does not handle multi-arm trials correctly.")
  } else {
    model_ume <- mtc.model(
      network,
      type = "ume",
      n.chain = n_chain,
      likelihood = "binom",
      link = "logit",
      linearModel = "random"
    )

    set.seed(seed)
    results_ume <- mtc.run(
      model_ume,
      n.adapt = n_adapt,
      n.iter = n_iter,
      thin = thin
    )
  }

  # Export model summaries and core diagnostics.
  write_text_output(summary(results_consistency), file.path("results", paste0("summary_", slug, ".txt")))
  if (!is.null(results_ume)) {
    write_text_output(summary(results_ume), file.path("results", paste0("summary_ume_", slug, ".txt")))
  } else {
    writeLines(
      c(
        "UME model not run.",
        "Reason: gemtc warns that UME does not handle multi-arm trials correctly for this network."
      ),
      file.path("results", paste0("summary_ume_", slug, ".txt"))
    )
  }

  relative_effects <- relative.effect(
    results_consistency,
    reference_treatment,
    preserve.extra = FALSE
  )

  save_forest_plot(
    relative_effects,
    file.path("figures", paste0("forest_", slug, ".png")),
    digits = 3
  )

  save_png(
    file.path("figures", paste0("trace_", slug, ".png")),
    function() plot(results_consistency)
  )

  save_png(
    file.path("figures", paste0("gelman_", slug, ".png")),
    function() gelman.plot(results_consistency)
  )

  # Export ranking and league-table results.
  ranks <- rank.probability(results_consistency, preferredDirection = preferred_direction)

  save_png(
    file.path("figures", paste0("rank_probability_", slug, ".png")),
    function() plot(ranks, beside = TRUE)
  )

  save_cumulative_rank_plot(ranks, outcome_label, file.path("figures", paste0("sucra_", slug, ".png")))

  sucra_values <- sucra(ranks)
  sucra_results <- data.frame(
    treatment = names(sucra_values),
    sucra = as.numeric(sucra_values),
    row.names = NULL
  )

  write.csv(sucra_results, file.path("results", paste0("sucra_", slug, ".csv")), row.names = FALSE)

  league_table_log_or <- relative.effect.table(results_consistency)
  league_table_or <- round(exp(league_table_log_or), 2)

  write.csv(league_table_log_or, file.path("results", paste0("league_table_log_or_", slug, ".csv")))
  write.csv(league_table_or, file.path("results", paste0("league_table_or_", slug, ".csv")))
}


# Step 7. Run both example analyses from the script directory.
project_dir <- get_script_dir()
setwd(project_dir)

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

message("Working directory: ", project_dir)

run_continuous_nma(file.path(project_dir, "ADAS_Cog.csv"), "ADAS_Cog")
run_binary_nma(file.path(project_dir, "AEs.csv"), "AEs")

message("Finished. Outputs written to ", file.path(project_dir, "results"), " and ", file.path(project_dir, "figures"))