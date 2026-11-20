wrapper_fixture <- function() {
  sim <- simulate_gradient_regions(240L, minority = 0.25, seed = 246L)
  n <- nrow(sim$xy)
  ids <- paste0("cell", seq_len(2L * n))
  xy <- rbind(sim$xy, sim$xy)
  rownames(xy) <- ids
  colnames(xy) <- c("x", "y")
  clusters <- factor(c(as.character(sim$labels),
                       rev(as.character(sim$labels))), levels = c(levels(sim$labels), "unused"))
  list(ids = ids, xy = xy, clusters = clusters,
       samples = rep(c("one", "two"), each = n),
       counts = matrix(1L, nrow = 2, ncol = length(ids),
                       dimnames = list(c("gene1", "gene2"), ids)))
}

plain_labels <- function(x) {
  structure(as.integer(x), levels = levels(x), class = "factor")
}

test_that("the coordinate API is exported and object input is validated", {
  expect_identical(getExportedValue("fibermargin", "FiberMargin"), FiberMargin)
  expect_named(formals(FiberMargin), c("xy", "labels", "samples", "workers", "transport"))
  expect_error(RunFiberMargin(list()), "must be a Seurat")
  expect_error(RunFiberMargin(list(), output.column = NA_character_), "output.column")
})

test_that("SpatialExperiment writes the same 3D specimen-wise labels as FiberMargin", {
  skip_if_not_installed("SpatialExperiment")
  skip_if_not_installed("SummarizedExperiment")
  x <- wrapper_fixture()
  xyz <- cbind(x$xy, z = ifelse(x$samples == "one", 0, x$xy[, 1L]^2))
  metadata <- data.frame(clusters = x$clusters, sample_id = x$samples,
                         row.names = x$ids)
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = x$counts), spatialCoords = xyz, colData = metadata
  )
  direct <- FiberMargin(xyz, x$clusters, x$samples, workers = 1L)
  expect_identical(unname(attr(direct, "dimensions_used")), c(2L, 3L))
  out <- RunFiberMargin(spe, cluster.column = "clusters", workers = 1L)
  expect_s4_class(out, "SpatialExperiment")
  expect_identical(out$FiberMargin_clusters, plain_labels(direct))
  short <- RunFiberMargin(spe, cluster.column = "clusters", workers = 1L,
                          transport = "short")
  expect_identical(short$FiberMargin_clusters, plain_labels(
    FiberMargin(xyz, x$clusters, x$samples, workers = 1L, transport = "short")
  ))
  expect_identical(out$clusters, spe$clusters)
  expect_identical(colnames(out), x$ids)
  expect_identical(SpatialExperiment::spatialCoords(out), SpatialExperiment::spatialCoords(spe))
  expect_identical(SummarizedExperiment::assay(out), SummarizedExperiment::assay(spe))
  expect_error(RunFiberMargin(spe), "Supply `cluster.column`")
  expect_error(RunFiberMargin(spe, "clusters", sample.column = "missing"), "missing")
  expect_error(RunFiberMargin(spe, "clusters", output.column = "clusters"), "already exists")
  expect_error(RunFiberMargin(spe, "clusters", images = "one"), "Unexpected arguments")
})

test_that("Seurat matches shuffled centroids and isolates overlapping coordinate systems", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  x <- wrapper_fixture()
  object <- suppressWarnings(SeuratObject::CreateSeuratObject(x$counts))
  object$clusters <- x$clusters
  SeuratObject::Idents(object) <- "clusters"
  for (sample in unique(x$samples)) {
    shuffled <- rev(which(x$samples == sample))
    object <- suppressWarnings(`[[<-`(object, sample, value = SeuratObject::CreateFOV(
      as.data.frame(x$xy[shuffled, ]), type = "centroids", assay = "RNA",
      key = paste0(sample, "_")
    )))
  }
  out <- RunFiberMargin(object, workers = 1L)
  original_idents <- SeuratObject::Idents(object)
  direct <- FiberMargin(x$xy, original_idents, x$samples, workers = 1L)
  expect_s4_class(out, "Seurat")
  expect_identical(out[[]]$FiberMargin_clusters, plain_labels(direct))
  expect_identical(SeuratObject::Idents(out), original_idents)
  expect_identical(out[[]]$clusters, x$clusters)
  one <- RunFiberMargin(object, "clusters", images = "one",
                        output.column = "one_refined", workers = 1L)
  expect_identical(one[[]]$one_refined[x$samples == "two"],
                   x$clusters[x$samples == "two"])
  expect_identical(levels(one[[]]$one_refined), levels(x$clusters))
  grouped <- object
  grouped$sample_id <- rep(c("p1", "p2"), length.out = length(x$ids))
  split <- RunFiberMargin(grouped, "clusters", workers = 1L)
  direct <- FiberMargin(x$xy, x$clusters, interaction(x$samples, grouped$sample_id),
                         workers = 1L)
  expect_identical(split[[]]$FiberMargin_clusters, plain_labels(direct))
  expect_error(RunFiberMargin(object, images = "missing"), "existing Seurat")
  object <- suppressWarnings(`[[<-`(object, "copy", value = object[["one"]]))
  expect_error(RunFiberMargin(object), "share cell IDs")
  expect_identical(
    RunFiberMargin(object, images = c("one", "two"), workers = 1L)[[]]$FiberMargin_clusters,
    out[[]]$FiberMargin_clusters
  )
})

test_that("Seurat Visium uses the full-resolution image coordinates", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("SeuratObject")
  x <- wrapper_fixture()
  counts <- x$counts[, x$samples == "one"]
  clusters <- x$clusters[x$samples == "one"]
  xy <- x$xy[x$samples == "one", ]
  object <- suppressWarnings(SeuratObject::CreateSeuratObject(counts, assay = "Spatial"))
  object$clusters <- clusters
  image <- methods::new("VisiumV1", image = array(0, c(2, 2, 3)),
    scale.factors = Seurat::scalefactors(lowres = 0.4),
    coordinates = data.frame(imagerow = xy[, "y"], imagecol = xy[, "x"]),
    spot.radius = 0.01, assay = "Spatial", key = "slice_"
  )
  object[["slice"]] <- image
  out <- RunFiberMargin(object, "clusters", workers = 1L)
  direct <- FiberMargin(xy, clusters, workers = 1L)
  expect_identical(out[[]]$FiberMargin_clusters, plain_labels(direct))

  if (methods::isClass("VisiumV2")) {
    fov <- SeuratObject::CreateFOV(as.data.frame(xy), type = "centroids",
                                  assay = "Spatial", key = "binned_")
    object[["slice"]] <- NULL
    object[["binned"]] <- methods::new("VisiumV2", fov,
      image = array(0, c(2, 2, 3)), scale.factors = Seurat::scalefactors(lowres = 0.4)
    )
    out <- RunFiberMargin(object, "clusters", workers = 1L)
    expect_identical(out[[]]$FiberMargin_clusters, plain_labels(direct))
  }
})

test_that("Giotto uses spatial location IDs, retains z, and writes new cell metadata", {
  skip_if_not_installed("GiottoClass")
  x <- wrapper_fixture()
  xyz <- cbind(x$xy, z = x$xy[, 1L]^2)
  loc <- data.frame(cell_ID = x$ids, sdimx = xyz[, 1L],
                    sdimy = xyz[, 2L], sdimz = xyz[, 3L])
  metadata <- data.frame(cell_ID = x$ids, clusters = x$clusters, sample_id = x$samples)
  object <- (function() {
    old <- options(giotto.use_conda = FALSE, giotto.no_python_warn = TRUE)
    on.exit(options(old))
    suppressWarnings(GiottoClass::createGiottoObject(
      expression = x$counts, spatial_locs = loc[rev(seq_len(nrow(loc))), ],
      cell_metadata = metadata, cores = 1L
    ))
  })()
  observed <- GiottoClass::getCellMetadata(object, output = "data.table")
  out <- RunFiberMargin(object, "clusters", workers = 1L)
  expect_s4_class(out, "giotto")
  refined <- GiottoClass::getCellMetadata(out, output = "data.table")
  direct <- FiberMargin(xyz, x$clusters, x$samples, workers = 1L)
  expect_identical(refined[["FiberMargin_clusters"]][match(x$ids, refined[["cell_ID"]])],
                   plain_labels(direct))
  expect_identical(refined[["clusters"]], observed[["clusters"]])
  expect_false("FiberMargin_clusters" %in%
                 colnames(GiottoClass::getCellMetadata(object, output = "data.table")))
})
