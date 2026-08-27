#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[[1L]] else
  file.path("benchmarks", "results", "release_0.1.0_worker_benchmark")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

library(fibermargin)

benchmark_data <- simulate_spatial_domains(
  n = 100000L,
  dimensions = 3L,
  pattern = "layers3d",
  noise = 0.20,
  noise_type = "boundary",
  samples = 4L,
  seed = 20260827L
)

run_once <- function(workers, repetition, measured = TRUE) {
  gc(verbose = FALSE)
  elapsed <- system.time({
    refined <- refine_spatial_labels(
      benchmark_data$xy,
      benchmark_data$labels,
      benchmark_data$samples,
      workers = workers
    )
  })[["elapsed"]]
  stopifnot(
    length(refined) == nrow(benchmark_data$xy),
    all(attr(refined, "dimensions_used") == 3L),
    sum(attr(refined, "sample_sizes")) == nrow(benchmark_data$xy)
  )
  data.frame(
    workers_requested = workers,
    workers_used = attr(refined, "workers"),
    repetition = repetition,
    measured = measured,
    elapsed_seconds = unname(elapsed),
    observations = nrow(benchmark_data$xy),
    specimens = nlevels(benchmark_data$samples),
    dimensions = ncol(benchmark_data$xy),
    stringsAsFactors = FALSE
  )
}

# Warm both native paths before timing.
invisible(run_once(1L, 0L, measured = FALSE))
invisible(run_once(4L, 0L, measured = FALSE))

schedule <- c(1L, 4L, 4L, 1L, 1L, 4L)
records <- do.call(rbind, lapply(seq_along(schedule), function(index) {
  run_once(schedule[[index]], index, measured = TRUE)
}))
summary <- aggregate(
  elapsed_seconds ~ workers_requested + workers_used + observations +
    specimens + dimensions,
  data = records,
  FUN = median
)
summary$speedup_vs_one <- summary$elapsed_seconds[summary$workers_requested == 1L] /
  summary$elapsed_seconds

write.csv(
  records,
  file.path(output_dir, "worker_timing_raw.csv"),
  row.names = FALSE
)
write.csv(
  summary,
  file.path(output_dir, "worker_timing_summary.csv"),
  row.names = FALSE
)

hardware <- c(
  sprintf("package_version: %s", as.character(packageVersion("fibermargin"))),
  sprintf("R_version: %s", R.version.string),
  sprintf("system: %s", paste(Sys.info()[c("sysname", "release", "machine")], collapse = " ")),
  sprintf("cpu: %s", Sys.info()[["machine"]]),
  sprintf("physical_cores_detected: %s", parallel::detectCores(logical = FALSE)),
  sprintf("logical_cores_detected: %s", parallel::detectCores(logical = TRUE)),
  "preprocessing_included: validation, specimen indexing, axis reduction, normalization, and all native refinement",
  "simulation_generation_included: no"
)
writeLines(hardware, file.path(output_dir, "hardware.txt"))

print(summary)
