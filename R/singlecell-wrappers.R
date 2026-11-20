#' Refine clustering in a spatial transcriptomics object
#'
#' Reads spatial coordinates and existing cluster assignments, calls
#' [RefineMargin()], and writes a new factor column to observation metadata.
#' The object class, original labels, and active Seurat identities are retained.
#'
#' @param object A Seurat, SpatialExperiment, or Giotto Suite object.
#' @param cluster.column Name of the input metadata column. The default
#'   `NULL` uses active identities for Seurat. Supply a column for other objects.
#' @param output.column Name of a new metadata column for the refined labels.
#' @param sample.column Optional metadata column identifying independent
#'   specimens. The default uses `sample_id` when present. Seurat images
#'   are always independent; the specimen column can split them further.
#' @param workers Optional total CPU budget passed to [RefineMargin()].
#' @param ... Object-specific arguments. For Seurat, `images` selects image or
#'   field-of-view names. By default all images are used, and they must cover
#'   every observation exactly once. With an explicit selection, unselected
#'   observations keep their input labels. For Giotto, `spat_unit`,
#'   `feat_type`, and `location` select the spatial unit, feature type,
#'   and spatial-location layer; `NULL` uses Giotto's defaults.
#' @param transport Multiclass transport preset passed to [RefineMargin()].
#' @return The updated object. Use [RefineMargin()] directly to obtain
#'   pointwise and sample-level diagnostics.
#' @export
RunFiberMargin <- function(object, cluster.column = NULL,
                           output.column = "FiberMargin_clusters",
                           sample.column = NULL, workers = NULL, ...,
                           transport = "standard") {
  .fm_column_name(cluster.column, "cluster.column", optional = TRUE)
  .fm_column_name(sample.column, "sample.column", optional = TRUE)
  .fm_column_name(output.column, "output.column")
  UseMethod("RunFiberMargin")
}

#' @export
RunFiberMargin.default <- function(object, cluster.column = NULL,
                                    output.column = "FiberMargin_clusters",
                                    sample.column = NULL, workers = NULL, ...,
                                    transport = "standard") {
  stop("`object` must be a Seurat, SpatialExperiment, or Giotto Suite object.",
       call. = FALSE)
}

#' @export
RunFiberMargin.SpatialExperiment <- function(object, cluster.column = NULL,
                                             output.column = "FiberMargin_clusters",
                                             sample.column = NULL, workers = NULL, ...,
                                             transport = "standard") {
  if (...length()) stop("Unexpected arguments for SpatialExperiment.", call. = FALSE)
  .fm_require("SpatialExperiment")
  .fm_require("SummarizedExperiment")
  metadata <- SummarizedExperiment::colData(object)
  .fm_new_column(output.column, colnames(metadata))
  clusters <- .fm_clusters(metadata, cluster.column)
  samples <- .fm_samples(metadata, sample.column, default = "sample_id")
  # spatialCoords() is column-aligned by the SpatialExperiment class.
  xy <- SpatialExperiment::spatialCoords(object)
  metadata[[output.column]] <- .fm_refined_factor(xy, clusters, samples, workers, transport)
  SummarizedExperiment::colData(object) <- metadata
  object
}

#' @export
RunFiberMargin.Seurat <- function(object, cluster.column = NULL,
                                  output.column = "FiberMargin_clusters",
                                  sample.column = NULL, workers = NULL, ...,
                                  transport = "standard", images = NULL) {
  if (...length()) stop("Unexpected arguments for Seurat; use named `images`.",
                        call. = FALSE)
  .fm_require("Seurat")
  .fm_require("SeuratObject")
  ids <- colnames(object)
  metadata <- object[[]]
  .fm_new_column(output.column, c(colnames(metadata), names(object)))
  available <- SeuratObject::Images(object)
  if (!length(available)) stop("The Seurat object has no spatial images or fields of view.",
                               call. = FALSE)
  if (!is.null(images) &&
      (!is.character(images) || !length(images) || anyNA(images) ||
       anyDuplicated(images) || !all(images %in% available))) {
    stop("`images` must contain unique names of existing Seurat images or fields of view.",
         call. = FALSE)
  }
  selected <- if (is.null(images)) available else images
  coordinates <- lapply(selected, function(name) {
    .fm_seurat_coordinates(object[[name]], name)
  })
  observed <- unlist(lapply(coordinates, rownames), use.names = FALSE)
  rows <- match(observed, ids)
  if (anyNA(rows)) stop("Seurat spatial coordinates include unknown cell IDs.",
                        call. = FALSE)
  if (anyDuplicated(rows)) {
    stop("Selected Seurat images share cell IDs. Choose non-overlapping `images`.",
         call. = FALSE)
  }
  if (is.null(images) && length(rows) != length(ids)) {
    stop("Not every Seurat observation has coordinates. Use `images` to select a subset.",
         call. = FALSE)
  }
  if (any(vapply(coordinates, ncol, integer(1L)) == 3L)) {
    coordinates <- lapply(coordinates, function(x) {
      if (ncol(x) == 2L) cbind(x, z = 0) else x
    })
  }
  xy <- do.call(rbind, coordinates)
  metadata <- metadata[match(ids, rownames(metadata)), , drop = FALSE]
  if (is.null(cluster.column)) {
    identities <- SeuratObject::Idents(object)
    clusters <- identities[match(ids, names(identities))]
  } else {
    clusters <- .fm_clusters(metadata, cluster.column)
  }
  samples <- rep(seq_along(selected), vapply(coordinates, nrow, integer(1L)))
  subgroups <- .fm_samples(metadata, sample.column, default = "sample_id")
  if (!is.null(subgroups)) {
    subgroups <- subgroups[rows]
    if (.fiber_has_missing(subgroups)) {
      stop("`sample.column` must not contain missing identifiers for selected cells.",
           call. = FALSE)
    }
    samples <- interaction(samples, factor(subgroups), drop = TRUE, sep = ":")
  }
  clusters[rows] <- .fm_refined_factor(xy, clusters[rows], samples, workers, transport)
  names(clusters) <- ids
  SeuratObject::AddMetaData(object, metadata = clusters, col.name = output.column)
}

#' @export
RunFiberMargin.giotto <- function(object, cluster.column = NULL,
                                  output.column = "FiberMargin_clusters",
                                  sample.column = NULL, workers = NULL, ...,
                                  transport = "standard", spat_unit = NULL,
                                  feat_type = NULL, location = NULL) {
  if (...length()) stop("Unexpected arguments for Giotto; use named layer selectors.",
                        call. = FALSE)
  .fm_require("GiottoClass")
  metadata <- GiottoClass::getCellMetadata(object, spat_unit = spat_unit,
    feat_type = feat_type, output = "data.table", copy_obj = TRUE)
  .fm_new_column(output.column, colnames(metadata))
  ids <- .fm_ids(metadata[["cell_ID"]], nrow(metadata), "Giotto cell metadata")
  clusters <- .fm_clusters(metadata, cluster.column)
  samples <- .fm_samples(metadata, sample.column, default = "sample_id")
  coords <- GiottoClass::getSpatialLocations(object, spat_unit = spat_unit,
    name = location, output = "data.table", copy_obj = TRUE, verbose = FALSE)
  xy <- .fm_coordinate_matrix(coords, list(c("sdimx", "sdimy", "sdimz")),
                              "Giotto spatial locations")
  location_ids <- .fm_ids(coords[["cell_ID"]], nrow(xy), "Giotto spatial locations")
  rows <- match(ids, location_ids)
  if (anyNA(rows)) stop("Giotto spatial locations do not cover every selected cell ID.",
                        call. = FALSE)
  refined <- .fm_refined_factor(xy[rows, , drop = FALSE], clusters, samples, workers, transport)
  added <- data.frame(cell_ID = ids, stringsAsFactors = FALSE)
  added[[output.column]] <- refined
  GiottoClass::addCellMetadata(object, spat_unit = spat_unit, feat_type = feat_type,
    new_metadata = added, by_column = TRUE, column_cell_ID = "cell_ID")
}

.fm_seurat_coordinates <- function(image, name) {
  source <- paste0("Seurat image '", name, "'")
  if (inherits(image, "FOV")) {
    boundaries <- SeuratObject::Boundaries(image)
    centers <- boundaries[vapply(boundaries, function(x) {
      inherits(image[[x]], "Centroids")
    }, logical(1L))]
    if (!length(centers)) stop(source, " does not contain a centroid boundary.",
                               call. = FALSE)
    if (length(centers) > 1L && !"centroids" %in% centers) {
      stop(source, " has multiple centroid boundaries and none named 'centroids'.",
           call. = FALSE)
    }
    which <- if ("centroids" %in% centers) "centroids" else centers[[1L]]
    coords <- SeuratObject::GetTissueCoordinates(image, which = which, full = FALSE)
  } else {
    coords <- SeuratObject::GetTissueCoordinates(image, scale = NULL)
  }
  xy <- .fm_coordinate_matrix(coords,
    list(c("x", "y", "z"), c("imagecol", "imagerow")), source)
  ids <- if ("cell" %in% colnames(coords)) coords[["cell"]] else rownames(coords)
  rownames(xy) <- .fm_ids(ids, nrow(xy), source)
  xy
}

.fm_coordinate_matrix <- function(coords, axes, source) {
  coords <- as.data.frame(coords)
  for (allowed in axes) {
    if (all(allowed[1:2] %in% names(coords))) {
      xy <- as.matrix(coords[, intersect(allowed, names(coords)), drop = FALSE])
      if (is.numeric(xy)) return(xy)
    }
  }
  stop(source, " does not contain supported numeric spatial coordinate columns.",
       call. = FALSE)
}

.fm_ids <- function(ids, n, source) {
  ids <- as.character(ids)
  if (n == 0L || length(ids) != n || anyNA(ids) || any(!nzchar(ids)) ||
      anyDuplicated(ids)) {
    stop(source, " must provide one unique, non-missing cell ID per row.",
         call. = FALSE)
  }
  ids
}

.fm_require <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("This object requires the optional package '", package, "'.", call. = FALSE)
  }
}

.fm_column_name <- function(name, argument, optional = FALSE) {
  if (is.null(name) && optional) return(invisible(NULL))
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`", argument, "` must be one non-empty column name.", call. = FALSE)
  }
}

.fm_new_column <- function(name, existing) {
  if (name %in% existing) {
    stop("`output.column` already exists: '", name, "'. Supply a new name.",
         call. = FALSE)
  }
}

.fm_clusters <- function(metadata, column) {
  if (is.null(column)) stop("Supply `cluster.column` for this object.", call. = FALSE)
  as.factor(.fm_metadata_column(metadata, column))
}

.fm_samples <- function(metadata, column, default = NULL) {
  if (is.null(column) && !is.null(default) && default %in% colnames(metadata)) {
    column <- default
  }
  if (!is.null(column)) .fm_metadata_column(metadata, column)
}

.fm_metadata_column <- function(metadata, column) {
  if (!column %in% colnames(metadata)) {
    stop("Observation metadata does not contain column '", column, "'.",
         call. = FALSE)
  }
  values <- metadata[[column]]
  if (isS4(values)) values <- as.vector(values)
  if (!is.atomic(values) || !is.null(dim(values))) {
    stop("Observation metadata column '", column, "' must be a vector.",
         call. = FALSE)
  }
  values
}

.fm_refined_factor <- function(xy, clusters, samples, workers, transport) {
  refined <- RefineMargin(xy, clusters, samples = samples, workers = workers,
                          transport = transport)
  # Container metadata stores assignments, not per-point engine diagnostics.
  structure(as.integer(refined), levels = levels(refined), class = "factor")
}
