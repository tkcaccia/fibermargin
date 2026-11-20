#!/usr/bin/env Rscript

# Rscript benchmarks/benchmark_public_api.R [output-directory]
suppressPackageStartupMessages(library(fibermargin))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) <= 1L)
output <- if (length(args)) args[[1L]] else "fibermargin-api-benchmark"
dir.create(output, showWarnings = FALSE, recursive = TRUE)

cases <- list(
  jagged_2d = simulate_spatial_domains(
    n = 6000L, pattern = "jagged_stripes", noise_type = "boundary",
    noise = 0.25, density_profile = "extreme", seed = 810001L
  ),
  gradient_2d = simulate_gradient_regions(
    n = 6000L, minority = 0.25, density_profile = "strong", seed = 810002L
  ),
  volume_3d = simulate_volumetric_domains(
    n = 60000L, shape = "folded_layers", acquisition = "irregular_z", samples = 3L,
    noise_type = "region", seed = 810003L
  ),
  DLPFC = load_spatial_benchmark("dlpfc", scenario = 1L, seed = 810004L),
  CRC = load_spatial_benchmark("crc", scenario = "CRC_random_25_r1", seed = 1040043L)
)

# Timing in benchmark_spatial_refiners() covers validation and native work only.
# A separate replay checks exact label and diagnostic equality across worker counts.
for (name in names(cases)) {
  case <- cases[[name]]
  serial <- RefineMargin(case$xy, case$labels, case$samples, workers = 1L)
  parallel <- RefineMargin(case$xy, case$labels, case$samples, workers = 4L)
  stopifnot(identical(serial, refine_spatial_labels(
    case$xy, case$labels, case$samples, workers = 1L
  ), num.eq = FALSE))
  attr(serial, "workers") <- attr(parallel, "workers") <- NULL
  stopifnot(identical(serial, parallel, num.eq = FALSE))
  message(name, ": labels and diagnostics identical across workers and aliases")
}
metrics <- benchmark_spatial_refiners(cases, seed = 820001L)
stopifnot(nrow(metrics) == 2L * length(cases), all(is.na(metrics$error)),
          identical(unique(metrics$method), c("Initial", "FiberMargin")))
utils::write.csv(metrics, file.path(output, "metrics.csv"), row.names = FALSE)
writeLines(c(
  paste("fibermargin", utils::packageVersion("fibermargin")),
  "Timing: one call with the default CPU budget; generation and scoring excluded.",
  capture.output(sessionInfo())
), file.path(output, "session-info.txt"))
print(metrics[, c("dataset", "method", "n", "accuracy", "ari", "damage_rate", "seconds")],
      row.names = FALSE)
