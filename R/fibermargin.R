#' Refine noisy categorical labels on a spatial domain
#'
#' `refine_spatial_labels()` applies FiberMargin to a coordinate-indexed label
#' field. The function is deterministic: no training phase is run and all route
#' constants are fixed internally. For multiclass data, each class receives a
#' two-sided local evidence score from short directional neighborhood sweeps; in
#' the binary case it applies a fixed nearest-neighbour ballot rule.
#'
#' The method can process multiple specimens at once. `samples` defines
#' independent coordinate systems: coordinates and labels are never shared
#' across specimens, even when their numeric coordinate ranges overlap. Constant
#' coordinate axes are removed separately within each specimen. Consequently, a
#' three-column specimen with constant `z` is processed identically to its
#' two-column `(x, y)` representation, while genuinely three-dimensional
#' specimens retain all three axes.
#'
#' @param xy Numeric matrix with two or three spatial coordinates per row.
#' @param labels Initial categorical assignment for every row of `xy`.
#' @param samples Optional specimen identifier. Different specimens are refined
#'   independently.
#' @param workers Optional total CPU budget shared by the native computation.
#'   `NULL` uses up to two physical cores. The budget is never multiplied by the
#'   number of specimens, and results are independent of the worker count.
#'
#' @return A factor with the levels and names of `labels`; when `labels` is
#'   unnamed, row names from `xy` are used as names. Attributes
#'   `candidate`, `margin_score`, `required`, `repair_margin`,
#'   `atlas_dispersion`, `isolation`, and `changed` contain pointwise
#'   diagnostics. `margin_score` is a local support contrast on the internal
#'   coordinate lattice, `repair_margin` is the difference between
#'   `margin_score` and the adaptive admission threshold, and `changed`
#'   marks updated sites. `isolation` is a deterministic local-gap
#'   protection factor in the multiclass route (one for the binary ballot). A
#'   candidate is accepted exactly when it differs from the observed label and
#'   `repair_margin` is nonnegative. Summary attributes report `workers`,
#'   `dimensions_used`, `labels_changed`, `changed_fraction`, `classes_before`,
#'   `classes_after`, `removed_classes`, and `sample_sizes`. Except for
#'   `workers`, these summaries are named by specimen; class summaries are named
#'   lists. Class counts are not constrained, so `removed_classes` explicitly
#'   records observed classes that disappear after refinement.
#' @export
#' @examples
#' sim <- simulate_gradient_regions(n = 2000, samples = 2)
#' refined <- refine_spatial_labels(sim$xy, sim$labels, sim$samples)
#' mean(refined == sim$truth)
#'
#' # A constant z axis is removed within this specimen.
#' xyz <- cbind(sim$xy, z = 0)
#' stopifnot(identical(
#'   refine_spatial_labels(xyz, sim$labels, sim$samples, workers = 1),
#'   refine_spatial_labels(sim$xy, sim$labels, sim$samples, workers = 1)
#' ))
refine_spatial_labels <- function(xy, labels, samples = NULL, workers = NULL) {
  .fiber_margin_engine(xy, labels, samples, workers = workers, control = list())
}

.fiber_margin_engine <- function(xy, labels, samples = NULL, workers = NULL,
                                 control = list()) {
  xy <- as.matrix(xy)
  if (!is.numeric(xy) || length(dim(xy)) != 2L) {
    stop("`xy` must be a numeric matrix.", call. = FALSE)
  }
  if (!ncol(xy) %in% c(2L, 3L)) {
    stop("`xy` must have exactly two or three coordinate columns.",
         call. = FALSE)
  }
  if (nrow(xy) < 1L) {
    stop("`xy` must contain at least one row.", call. = FALSE)
  }
  if (any(!is.finite(xy))) {
    stop("`xy` must contain only finite coordinates.", call. = FALSE)
  }
  storage.mode(xy) <- "double"
  if (length(labels) != nrow(xy)) {
    stop("`labels` must have one entry per row of `xy`.", call. = FALSE)
  }
  if (anyNA(labels)) {
    stop("`labels` must not contain missing assignments.", call. = FALSE)
  }
  output_names <- names(labels)
  if (is.null(output_names)) output_names <- rownames(xy)
  labels <- as.factor(labels)
  if (is.null(samples)) {
    samples <- rep.int(1L, nrow(xy))
  } else if (!(is.factor(samples) || is.character(samples) ||
               is.integer(samples) ||
               (is.numeric(samples) && all(is.finite(samples)) &&
                all(samples == floor(samples))))) {
    stop("`samples` must be `NULL`, integer, character, or factor.",
         call. = FALSE)
  }
  if (length(samples) != nrow(xy)) {
    stop("`samples` must have one identifier per row of `xy`.", call. = FALSE)
  }
  if (anyNA(samples)) {
    stop("`samples` must not contain missing identifiers.", call. = FALSE)
  }
  samples <- droplevels(as.factor(samples))
  if (!is.list(control) || (length(control) && is.null(names(control)))) {
    stop("Internal FiberMargin control must be a named list.", call. = FALSE)
  }

  if (is.null(workers)) {
    cores <- parallel::detectCores(logical = FALSE)
    if (is.na(cores)) cores <- 1L
    workers <- max(1L, min(2L, as.integer(cores) - 1L))
  }
  if (!is.numeric(workers) || length(workers) != 1L ||
      !is.finite(workers) || workers < 1 || workers != floor(workers) ||
      workers > .Machine$integer.max) {
    stop("`workers` must be `NULL` or one positive integer.", call. = FALSE)
  }
  workers <- as.integer(workers)
  if (tolower(Sys.getenv("_R_CHECK_LIMIT_CORES_", "false")) %in%
      c("true", "yes") && workers > 2L) {
    workers <- 2L
  }

  groups <- split(seq_len(nrow(xy)), samples, drop = TRUE)
  dimensions_used <- vapply(groups, function(rows) {
    sum(vapply(seq_len(ncol(xy)), function(axis) {
      values <- xy[rows, axis]
      min(values) < max(values)
    }, logical(1L)))
  }, integer(1L))
  insufficient <- which(dimensions_used < 2L)
  if (length(insufficient)) {
    sample_name <- names(dimensions_used)[insufficient[1L]]
    stop(sprintf(
      "Sample `%s` has fewer than two varying coordinate axes after constant axes are removed.",
      sample_name
    ), call. = FALSE)
  }

  native_control <- control
  native_control$threads <- workers
  result <- .Call(
    `_fibermargin_fiber_margin_cpp`,
    xy,
    as.integer(labels),
    as.integer(samples),
    native_control
  )
  native_sample_names <- levels(samples)[result$sample_codes]
  names(result$dimensions_used) <- native_sample_names
  names(result$sample_sizes) <- native_sample_names
  if (!identical(unname(result$dimensions_used),
                 unname(dimensions_used[native_sample_names]))) {
    stop("Internal dimensionality diagnostics are inconsistent.", call. = FALSE)
  }

  refined <- factor(levels(labels)[result$labels], levels = levels(labels))
  names(refined) <- output_names
  candidate <- factor(
    levels(labels)[result$candidate], levels = levels(labels)
  )
  names(candidate) <- output_names
  pointwise <- list(
    candidate = candidate,
    margin_score = result$margin_score,
    required = result$required,
    repair_margin = result$margin_score - result$required,
    atlas_dispersion = result$atlas_dispersion,
    isolation = result$isolation,
    changed = result$changed
  )
  for (name in names(pointwise)[-1L]) names(pointwise[[name]]) <- output_names
  for (name in names(pointwise)) attr(refined, name) <- pointwise[[name]]

  classes_for <- function(values) {
    levels(labels)[sort(unique(as.integer(values)))]
  }
  classes_before <- lapply(groups, function(rows) classes_for(labels[rows]))
  classes_after <- lapply(groups, function(rows) classes_for(refined[rows]))
  removed_classes <- Map(setdiff, classes_before, classes_after)
  labels_changed <- vapply(groups, function(rows) {
    sum(result$changed[rows])
  }, integer(1L))
  sample_sizes <- lengths(groups)
  changed_fraction <- labels_changed / sample_sizes

  attr(refined, "workers") <- workers
  attr(refined, "dimensions_used") <- result$dimensions_used
  attr(refined, "labels_changed") <- labels_changed
  attr(refined, "changed_fraction") <- changed_fraction
  attr(refined, "classes_before") <- classes_before
  attr(refined, "classes_after") <- classes_after
  attr(refined, "removed_classes") <- removed_classes
  attr(refined, "sample_sizes") <- result$sample_sizes
  refined
}
