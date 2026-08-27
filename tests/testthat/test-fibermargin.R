test_that("FiberMargin returns stable factors and diagnostics", {
  simulation <- simulate_spatial_domains(
    n = 1200, pattern = "jagged_stripes", noise = 0.20,
    noise_type = "boundary", seed = 1501
  )
  refined <- refine_spatial_labels(
    simulation$xy, simulation$labels, simulation$samples, workers = 1
  )

  expect_s3_class(refined, "factor")
  expect_identical(levels(refined), levels(factor(simulation$labels)))
  expect_length(refined, nrow(simulation$xy))
  expect_length(attr(refined, "candidate"), nrow(simulation$xy))
  expect_true(all(is.finite(attr(refined, "margin_score"))))
  expect_true(all(is.finite(attr(refined, "required"))))
  expect_equal(
    attr(refined, "repair_margin"),
    attr(refined, "margin_score") - attr(refined, "required")
  )
  expect_null(attr(refined, "risk"))
  expect_true(all(attr(refined, "atlas_dispersion") >= 0))
  expect_identical(
    unname(attr(refined, "changed")),
    unname(refined != factor(simulation$labels))
  )

  observed <- as.integer(factor(simulation$labels, levels = levels(refined)))
  candidate <- as.integer(attr(refined, "candidate"))
  expect_true(all(attr(refined, "required") >= 0))
  expect_identical(
    attr(refined, "changed"),
    candidate != observed & attr(refined, "repair_margin") >= 0
  )
})

test_that("binary local-volume evidence uses the shared margin admission rule", {
  set.seed(1551)
  xy <- matrix(runif(1600), ncol = 2)
  labels <- factor(ifelse(xy[, 1] + 0.12 * sin(8 * xy[, 2]) < 0.5, "A", "B"))
  noisy <- labels
  corrupted <- sample.int(length(noisy), 120L)
  noisy[corrupted] <- ifelse(noisy[corrupted] == "A", "B", "A")
  noisy <- factor(noisy, levels = levels(labels))

  refined <- refine_spatial_labels(xy, noisy, workers = 1)
  observed <- as.integer(noisy)
  candidate <- as.integer(attr(refined, "candidate"))

  expect_true(all(attr(refined, "required") == 1))
  expect_true(all(attr(refined, "margin_score")[candidate == observed] == 0))
  expect_true(all(candidate[attr(refined, "margin_score") == 0] ==
                  observed[attr(refined, "margin_score") == 0]))
  expect_identical(
    attr(refined, "changed"),
    candidate != observed & attr(refined, "repair_margin") >= 0
  )
})

test_that("overlapping specimens are exactly isolated", {
  set.seed(1601)
  shared_xy <- matrix(runif(2400), ncol = 2)
  first_labels <- cut(
    shared_xy[, 1] + 0.1 * sin(8 * shared_xy[, 2]),
    breaks = c(-Inf, 0.33, 0.66, Inf), labels = c("A", "B", "C")
  )
  second_labels <- cut(
    shared_xy[, 2], breaks = c(-Inf, 0.33, 0.66, Inf),
    labels = c("A", "B", "C")
  )
  xy <- rbind(shared_xy, shared_xy)
  labels <- factor(
    c(as.character(first_labels), as.character(second_labels)),
    levels = c("A", "B", "C", "D")
  )
  samples <- factor(rep(c("first", "second"), each = nrow(shared_xy)))

  original <- refine_spatial_labels(xy, labels, samples, workers = 4)
  modified_labels <- labels
  modified_labels[samples == "second"] <- sample(
    levels(labels), sum(samples == "second"), replace = TRUE
  )
  modified <- refine_spatial_labels(
    xy, modified_labels, samples, workers = 2
  )
  first <- samples == "first"

  expect_identical(
    as.character(original[first]), as.character(modified[first])
  )
  for (diagnostic in c(
    "candidate", "margin_score", "required", "repair_margin",
    "atlas_dispersion", "isolation", "changed"
  )) {
    expect_identical(
      unname(attr(original, diagnostic)[first]),
      unname(attr(modified, diagnostic)[first])
    )
  }
})

test_that("chart and specimen workers are bitwise deterministic", {
  simulation <- simulate_spatial_domains(
    n = 2400, pattern = "jagged_stripes", noise = 0.25,
    samples = 3, seed = 1651
  )
  serial <- refine_spatial_labels(
    simulation$xy, simulation$labels, simulation$samples, workers = 1
  )
  parallel_two <- refine_spatial_labels(
    simulation$xy, simulation$labels, simulation$samples, workers = 2
  )
  parallel_four <- refine_spatial_labels(
    simulation$xy, simulation$labels, simulation$samples, workers = 4
  )

  expect_identical(attr(serial, "workers"), 1L)
  expect_lte(attr(parallel_two, "workers"), 2L)
  expect_lte(attr(parallel_four, "workers"), 4L)
  attr(serial, "workers") <- NULL
  attr(parallel_two, "workers") <- NULL
  attr(parallel_four, "workers") <- NULL
  expect_identical(parallel_two, serial)
  expect_identical(parallel_four, serial)

  single <- simulate_spatial_domains(
    n = 1600, pattern = "rings", noise = 0.25, samples = 1, seed = 1652
  )
  serial_chart <- refine_spatial_labels(
    single$xy, single$labels, single$samples, workers = 1
  )
  parallel_chart <- refine_spatial_labels(
    single$xy, single$labels, single$samples, workers = 2
  )
  attr(serial_chart, "workers") <- NULL
  attr(parallel_chart, "workers") <- NULL
  expect_identical(parallel_chart, serial_chart)
})

test_that("CRAN core limits cap explicit worker budgets", {
  old_limit <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA_character_)
  on.exit({
    if (is.na(old_limit)) {
      Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
    } else {
      Sys.setenv("_R_CHECK_LIMIT_CORES_" = old_limit)
    }
  })
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "true")

  simulation <- simulate_spatial_domains(
    n = 600, pattern = "rings", noise = 0.20, samples = 2, seed = 1653
  )
  refined <- refine_spatial_labels(
    simulation$xy, simulation$labels, simulation$samples, workers = 8
  )

  expect_lte(attr(refined, "workers"), 2L)
})

test_that("samples without two varying axes are rejected", {
  xy <- matrix(0, nrow = 30, ncol = 2)
  labels <- factor(rep(c("A", "B", "C"), length.out = 30))

  expect_error(
    refine_spatial_labels(xy, labels, workers = 1),
    "fewer than two varying coordinate axes"
  )
})

test_that("constant axes are removed separately within each specimen", {
  simulation <- simulate_spatial_domains(
    n = 1800, pattern = "jagged_stripes", noise = 0.20,
    samples = 3, seed = 1661
  )
  two_dimensional <- refine_spatial_labels(
    simulation$xy, simulation$labels, simulation$samples, workers = 1
  )
  constant_z <- refine_spatial_labels(
    cbind(simulation$xy, z = 7), simulation$labels,
    simulation$samples, workers = 1
  )

  expect_identical(constant_z, two_dimensional)
  expect_true(all(attr(constant_z, "dimensions_used") == 2L))

  set.seed(1662)
  n <- 300L
  mixed_xy <- rbind(
    cbind(matrix(runif(2L * n), ncol = 2), z = 0),
    matrix(runif(3L * n), ncol = 3)
  )
  mixed_labels <- factor(rep(c("A", "B", "C"), length.out = 2L * n))
  mixed_samples <- rep(c("flat", "volume"), each = n)
  mixed <- refine_spatial_labels(
    mixed_xy, mixed_labels, mixed_samples, workers = 2
  )
  expect_identical(
    attr(mixed, "dimensions_used"), c(flat = 2L, volume = 3L)
  )
})

test_that("genuine three-dimensional structure is not flattened", {
  set.seed(1663)
  base_xy <- matrix(runif(800), ncol = 2)
  xyz <- cbind(
    base_xy[rep(seq_len(nrow(base_xy)), 3L), , drop = FALSE],
    z = rep(c(-1, 0, 1), each = nrow(base_xy))
  )
  truth <- factor(rep(c("A", "B", "C"), each = nrow(base_xy)))
  noisy <- truth
  corrupted <- sample.int(length(noisy), 180L)
  noisy[corrupted] <- sample(levels(noisy), length(corrupted), replace = TRUE)
  noisy <- factor(noisy, levels = levels(truth))

  refined_3d <- refine_spatial_labels(xyz, noisy, workers = 1)
  refined_2d <- refine_spatial_labels(xyz[, 1:2], noisy, workers = 1)

  expect_identical(unname(attr(refined_3d, "dimensions_used")), 3L)
  expect_gt(sum(refined_3d != refined_2d), 0L)
  expect_gt(mean(refined_3d == truth), mean(refined_2d == truth))
})

test_that("row permutation is equivariant", {
  simulation <- simulate_spatial_domains(
    n = 1800, pattern = "branching", noise = 0.22,
    samples = 3, seed = 1664
  )
  baseline <- refine_spatial_labels(
    simulation$xy, simulation$labels, simulation$samples, workers = 1
  )
  set.seed(1665)
  permutation <- sample(nrow(simulation$xy))
  inverse <- order(permutation)
  permuted <- refine_spatial_labels(
    simulation$xy[permutation, , drop = FALSE],
    simulation$labels[permutation], simulation$samples[permutation],
    workers = 2
  )

  expect_identical(
    as.character(permuted[inverse]), as.character(baseline)
  )
  for (diagnostic in c(
    "candidate", "margin_score", "required", "repair_margin",
    "atlas_dispersion", "isolation", "changed"
  )) {
    expect_identical(
      unname(attr(permuted, diagnostic)[inverse]),
      unname(attr(baseline, diagnostic))
    )
  }
})

test_that("translation and uniform scaling preserve the result", {
  simulation <- simulate_spatial_domains(
    n = 1500, pattern = "rings", noise = 0.20,
    samples = 2, seed = 1666
  )
  baseline <- refine_spatial_labels(
    simulation$xy, simulation$labels, simulation$samples, workers = 1
  )
  transformed <- refine_spatial_labels(
    17 + 3.5 * simulation$xy, simulation$labels,
    simulation$samples, workers = 2
  )
  attr(baseline, "workers") <- NULL
  attr(transformed, "workers") <- NULL
  expect_identical(transformed, baseline)
})

test_that("integer character and factor inputs agree", {
  simulation <- simulate_spatial_domains(
    n = 1200, pattern = "jagged_stripes", noise = 0.20,
    samples = 3, seed = 1667
  )
  label_codes <- as.integer(simulation$labels)
  sample_codes <- as.integer(factor(simulation$samples))
  integer_result <- refine_spatial_labels(
    simulation$xy, label_codes, sample_codes, workers = 1
  )
  character_result <- refine_spatial_labels(
    simulation$xy, as.character(label_codes), as.character(sample_codes),
    workers = 1
  )
  factor_result <- refine_spatial_labels(
    simulation$xy, factor(as.character(label_codes)),
    factor(as.character(sample_codes)), workers = 1
  )

  expect_identical(integer_result, character_result)
  expect_identical(integer_result, factor_result)
})

test_that("names levels and class-loss summaries are preserved", {
  xy <- as.matrix(expand.grid(x = seq_len(20L), y = seq_len(20L)))
  rownames(xy) <- paste0("cell_", seq_len(nrow(xy)))
  labels <- factor(rep("A", nrow(xy)), levels = c("unused", "A", "B"))
  labels[which.min(rowSums((xy - 10)^2))] <- "B"
  names(labels) <- rownames(xy)

  refined <- refine_spatial_labels(xy, labels, workers = 1)

  expect_identical(levels(refined), levels(labels))
  expect_identical(names(refined), names(labels))
  expect_identical(names(attr(refined, "candidate")), names(labels))
  expect_identical(attr(refined, "classes_before")[[1L]], c("A", "B"))
  expect_identical(attr(refined, "classes_after")[[1L]], "A")
  expect_identical(attr(refined, "removed_classes")[[1L]], "B")
  expect_false("unused" %in% attr(refined, "classes_before")[[1L]])
  expect_identical(attr(refined, "labels_changed"), c(`1` = 1L))
  expect_identical(attr(refined, "sample_sizes"), c(`1` = 400L))
})

test_that("input validation reports each contract violation", {
  xy <- matrix(runif(60), ncol = 2)
  labels <- factor(rep(c("A", "B", "C"), length.out = nrow(xy)))

  expect_error(refine_spatial_labels(xy[, 1, drop = FALSE], labels),
               "exactly two or three")
  expect_error(refine_spatial_labels(cbind(xy, xy), labels),
               "exactly two or three")
  nonfinite <- xy
  nonfinite[1, 1] <- Inf
  expect_error(refine_spatial_labels(nonfinite, labels), "finite coordinates")
  expect_error(refine_spatial_labels(xy, labels[-1]), "one entry per row")
  missing_labels <- labels
  missing_labels[1] <- NA
  expect_error(refine_spatial_labels(xy, missing_labels), "missing assignments")
  expect_error(refine_spatial_labels(xy, labels, samples = 1:2),
               "one identifier per row")
  missing_samples <- rep("sample", nrow(xy))
  missing_samples[1] <- NA
  expect_error(refine_spatial_labels(xy, labels, missing_samples),
               "missing identifiers")
  for (invalid in list(0, -1, 1.5, Inf, NA_real_, "2", c(1, 2))) {
    expect_error(refine_spatial_labels(xy, labels, workers = invalid),
                 "positive integer")
  }
  line_xy <- cbind(x = seq_len(nrow(xy)), y = 0)
  expect_error(refine_spatial_labels(line_xy, labels),
               "fewer than two varying")
})

test_that("large multi-sample execution respects one CPU budget", {
  simulation <- simulate_spatial_domains(
    n = 12000, pattern = "jagged_stripes", noise = 0.20,
    samples = 8, seed = 1668
  )
  refined <- refine_spatial_labels(
    simulation$xy, simulation$labels, simulation$samples, workers = 4
  )

  expect_length(refined, 12000L)
  expect_lte(attr(refined, "workers"), 4L)
  expect_identical(sum(attr(refined, "sample_sizes")), 12000L)
  expect_identical(length(attr(refined, "sample_sizes")), 8L)
  expect_lt(as.numeric(object.size(refined)), 10 * 1024^2)
})

test_that("an unused factor level is not introduced in a regular spatial sample", {
  set.seed(1701)
  xy <- matrix(runif(800), ncol = 2)
  labels <- factor(ifelse(xy[, 1] < 0.5, "A", "B"), levels = c("C", "A", "B"))
  refined <- refine_spatial_labels(xy, labels, workers = 1)
  expect_false(any(refined == "C"))
  expect_false(any(attr(refined, "candidate") == "C"))
})

test_that("FiberMargin exposes no mathematical tuning controls", {
  set.seed(1801)
  xy <- matrix(runif(82), ncol = 2)
  labels <- factor(c(rep("A", 20), rep("B", 21)), levels = c("A", "B"))
  expect_error(
    fibermargin:::.fiber_margin_engine(
      xy, labels, workers = 1,
      control = list(tail_fraction = 0.50)
    ),
    "no tuning controls"
  )
})

test_that("fixed-k modal reference uses the shared neighbour-mode kernel", {
  xy <- as.matrix(expand.grid(x = seq_len(5L), y = seq_len(5L)))
  labels <- factor(ifelse(xy[, 1L] <= 3L, "A", "B"))
  samples <- factor(rep("mask", nrow(xy)))

  modal <- fibermargin:::.local_modal_filter_labels(
    xy, labels, samples, neighbors = 8L
  )
  direct <- fibermargin:::.refine_published_labels(
    xy, labels, samples, method = "graphst", neighbors = 8L
  )

  expect_identical(modal, direct)
  expect_identical(attr(modal, "neighbors"), attr(direct, "neighbors"))
})

test_that("alpha-expansion Potts control preserves factor labels and diagnostics", {
  simulation <- simulate_spatial_domains(
    n = 800, pattern = "rings", noise = 0.20,
    noise_type = "random", seed = 1851
  )
  refined <- fibermargin:::.alpha_expansion_potts_labels(
    simulation$xy, simulation$labels, simulation$samples
  )

  expect_s3_class(refined, "factor")
  expect_identical(levels(refined), levels(simulation$labels))
  expect_length(refined, nrow(simulation$xy))
  expect_true(all(!is.na(refined)))
  expect_true(is.finite(attr(refined, "energy")))
  expect_true(attr(refined, "energy") >= 0)
  expect_true(all(attr(refined, "neighbors") >= 0L))
  expect_true(all(attr(refined, "changes") >= 0L))
})
