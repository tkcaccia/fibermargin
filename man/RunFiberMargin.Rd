% Please edit documentation in R/singlecell-wrappers.R
\name{RunFiberMargin}
\alias{RunFiberMargin}
\title{Refine clustering in a spatial transcriptomics object}
\usage{
RunFiberMargin(
  object,
  cluster.column = NULL,
  output.column = "FiberMargin_clusters",
  sample.column = NULL,
  workers = NULL,
  ...
)
}
\arguments{
\item{object}{A Seurat, SpatialExperiment, or Giotto Suite object.}
\item{cluster.column}{Name of the input metadata column. The default
\code{NULL} uses active identities for Seurat. Supply a column for other objects.}
\item{output.column}{Name of a new metadata column for the refined labels.}
\item{sample.column}{Optional column identifying independent specimens. The
default uses \code{sample_id} when present. Seurat images are always
independent; the specimen column can split them further.}
\item{workers}{Optional total CPU budget passed to \code{\link{RefineMargin}()}.}
\item{...}{Object-specific arguments. For Seurat, \code{images} selects image
or field-of-view names. By default all images must cover every observation
exactly once. With an explicit selection, unselected observations keep their
input labels. For Giotto, \code{spat_unit}, \code{feat_type}, and \code{location}
select the spatial unit, feature type and spatial-location layer; \code{NULL}
uses Giotto's defaults.}
}
\description{
Reads observation coordinates and existing labels, calls
\code{\link{RefineMargin}()}, and saves a factor in a new metadata column.
Input metadata, object class, and active Seurat identities are not changed.
}
\details{
Seurat uses spot coordinates or field-of-view centroids, never polygons or
expression embeddings. SpatialExperiment uses \code{spatialCoords()}.
Giotto uses \code{sdimx}, \code{sdimy} and, when present, \code{sdimz} from
the selected location layer. Seurat and Giotto locations are matched by
observation ID. All selected observations must have valid coordinates and
cluster assignments.

The object adapters use optional packages: Seurat and SeuratObject,
SpatialExperiment and SummarizedExperiment, or GiottoClass (from Giotto Suite).
The matrix API does not require them.
}
\value{
The updated object with the new column \code{output.column}.
Use \code{\link{RefineMargin}()} directly for pointwise and specimen-level
diagnostics.
}
\examples{
if (requireNamespace("SpatialExperiment", quietly = TRUE)) {
  sim <- simulate_gradient_regions(n = 300, samples = 2L, seed = 1L)
  spots <- paste0("spot", seq_len(nrow(sim$xy)))
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = matrix(1L, nrow = 1L, ncol = length(spots),
                                 dimnames = list("gene", spots))),
    spatialCoords = sim$xy,
    colData = data.frame(clusters = sim$labels, sample_id = sim$samples,
                         row.names = spots)
  )
  spe <- RunFiberMargin(spe, cluster.column = "clusters", workers = 1L)
  table(spe$clusters, spe$FiberMargin_clusters)
}
}
