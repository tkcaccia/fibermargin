#!/usr/bin/env Rscript

# Use separate R processes/library paths for the reference and candidate.
# Rscript benchmarks/benchmark_native_performance.R prepare INPUT_DIR
# Rscript benchmarks/benchmark_native_performance.R run INPUT_DIR CASE WORKERS REPS OUTPUT_RDS
# Rscript benchmarks/benchmark_native_performance.R compare OLD_RDS NEW_RDS
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L)
mode <- match.arg(args[1L], c("prepare", "run", "compare", "geometry", "keys"))
folder <- args[2L]
if (mode %in% c("prepare", "run", "geometry")) {
  suppressPackageStartupMessages(library(fibermargin))
  # Allow the same runner to audit releases predating the public rename.
  api <- intersect(c("RefineMargin", "refine_spatial_labels"),
                   getNamespaceExports("fibermargin"))
  refine <- getExportedValue("fibermargin", api[[1L]])
}

if (mode == "prepare") {
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  set.seed(110241L)
  count <- 60000L
  xy <- cbind(runif(count), runif(count))
  labels <- factor(ifelse(xy[, 1L] < .5 + .1 * sin(9 * xy[, 2L]), "A", "B"))
  flip <- sample.int(count, count %/% 4L)
  labels[flip] <- levels(labels)[3L - as.integer(labels[flip])]
  crc <- load_spatial_benchmark("crc", "CRC_random_25_r1", seed = 1040043L)
  make_volume <- function(n, k, specimens) {
    xy <- matrix(runif(3L * n, -1, 1), n, 3L)
    latent <- if (specimens == 1L) sqrt(rowSums(xy^2)) else
      xy[, 3L] + .2 * sin(6 * xy[, 1L]) * cos(6 * xy[, 2L])
    truth <- as.integer(cut(latent, quantile(latent, seq(0, 1, length.out = k + 1L)),
                            include.lowest = TRUE))
    noisy <- truth
    flip <- sample.int(n, n %/% 4L)
    rival <- sample.int(k - 1L, length(flip), replace = TRUE)
    noisy[flip] <- rival + (rival >= truth[flip])
    list(xy = xy, labels = factor(noisy, levels = seq_len(k)),
         samples = rep(seq_len(specimens), length.out = n))
  }
  volume <- make_volume(100000L, 5L, 4L)
  large <- make_volume(250000L, 16L, 1L)
  batch <- list(xy = cbind(runif(90000L), runif(90000L)),
                labels = factor(rep(1:3, 30000L)), samples = rep(1:6000, each = 15L))
  cases <- list(binary = list(xy = xy, labels = labels), crc = crc, volume = volume,
                large = large, batch = batch)
  for (name in names(cases)) saveRDS(cases[[name]], file.path(folder, paste0(name, ".rds")))
} else if (mode == "compare") {
  stopifnot(length(args) == 3L)
  old <- readRDS(folder)
  new <- readRDS(args[3L])
  stopifnot(identical(old, new, num.eq = FALSE))
  cat("Labels and every attached diagnostic match bit for bit.\n")
} else if (mode == "keys") {
  # For this mode supply OLD_CPP and NEW_CPP instead of INPUT_DIR.
  stopifnot(length(args) == 3L)
  extract <- function(file, start) {
    lines <- readLines(file)
    first <- grep(start, lines)
    last <- grep("^struct Geometry", lines)
    stopifnot(length(first) == 1L, length(last) == 1L, first < last)
    lines[first:(last - 1L)]
  }
  src <- c("// [[Rcpp::plugins(cpp17)]]", "#include <Rcpp.h>", "#include <algorithm>",
    "#include <array>", "#include <cstdint>", "#include <random>",
    "namespace original {", extract(folder, "^std::uint64_t hilbert_code"), "}",
    "namespace revised {", extract(args[3L], "^// Spread the 16 low bits"), "}",
    "// [[Rcpp::export]]",
    "void compare_keys(int count) {",
    "  std::mt19937_64 rng(774193);",
    "  for (int i = 0; i < count; ++i) {",
    "    double point[3];",
    "    for (double& x : point) {",
    "      x = i % 2 == 0 ? (static_cast<double>(rng() & 65535) + 0.5) / 65535.0 :",
    "         std::generate_canonical<double, 53>(rng) * 1.02 - 0.01;",
    "    }",
    "    if (original::hilbert_code(point, 2) != revised::hilbert_code<2>(point) ||",
    "        original::hilbert_code(point, 3) != revised::hilbert_code<3>(point))",
    "      Rcpp::stop(\"Hilbert key mismatch at point %d\", i);",
    "  }",
    "}")
  Rcpp::sourceCpp(code = paste(src, collapse = "\n"))
  compare_keys(1000000L)
  cat("One million points match in 2D and 3D.\n")
} else if (mode == "geometry") {
  # The PRL geometry design includes varying density, three corruption modes,
  # multi-specimen inputs, and both 2D and 3D cases.
  stopifnot(length(args) == 4L)
  design <- read.csv(args[3L])
  stopifnot(nrow(design) == 108L)
  results <- lapply(seq_len(nrow(design)), function(i) {
    x <- design[i, ]
    d <- simulate_complex_spatial_domains(5000L, x$shape, x$density_profile,
      x$noise_type, samples = x$samples, seed = 3100000L + x$scenario)
    refine(d$xy, d$labels, d$samples, workers = 1L)
  })
  saveRDS(results, args[4L])
} else {
  stopifnot(length(args) == 6L)
  case <- args[3L]
  workers <- as.integer(args[4L])
  repeats <- as.integer(args[5L])
  stopifnot(repeats > 0L)
  data <- readRDS(file.path(folder, paste0(case, ".rds")))
  run <- function() refine(data$xy, data$labels, data$samples, workers)
  gc(verbose = FALSE)
  result <- run()
  seconds <- vapply(seq_len(repeats), function(i) {
    system.time(invisible(run()))[["elapsed"]]
  }, numeric(1L))
  saveRDS(result, args[6L])
  write.csv(data.frame(case = case, workers = workers, seconds = seconds),
            stdout(), row.names = FALSE)
}
