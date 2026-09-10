# FiberMargin

[![R CMD check](https://github.com/tkcaccia/fibermargin/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tkcaccia/fibermargin/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/tkcaccia/fibermargin/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/tkcaccia/fibermargin/actions/workflows/pkgdown.yaml)

FiberMargin performs training-free repair of categorical masks and spatial labels.
It is a post-processing method: upstream segmentation, clustering, or annotation
remains unchanged. A deterministic atlas of rotated Hilbert paths transports
leave-self-out class evidence from both directions using normalized path length
only. Its class score is the two-sided enclosure
`(sqrt(left) + sqrt(right))^2`. FiberMargin changes a label only when its
full-atlas margin clears a conservative chart-disagreement barrier; locally
isolated points receive extra protection at that decision.

The package exposes one production refiner with four arguments. It does not need
expression values, a graph, training labels, or user-selected tuning parameters.

## Installation

```r
install.packages("remotes")
remotes::install_github("tkcaccia/fibermargin", build_vignettes = TRUE)
```

## Basic Use

```r
library(fibermargin)

sim <- simulate_gradient_regions(
  n = 50000,
  minority = 0.25,
  samples = 2,
  density_profile = "strong"
)

refined <- refine_spatial_labels(
  xy = sim$xy,
  labels = sim$labels,
  samples = sim$samples
)

c(
  input = mean(sim$labels == sim$truth),
  refined = mean(refined == sim$truth)
)
```

Only `xy` and `labels` are required. `samples` defines independent coordinate
systems: coordinates may overlap completely between specimens, but neither
coordinates nor labels cross a specimen boundary. Constant coordinate axes are
removed separately within each specimen. A specimen with columns `(x, y, z)`
and constant `z` is therefore exactly equivalent to `(x, y)`, while genuinely
three-dimensional specimens use all three coordinates.

`workers` is the total native CPU budget for the call, not a budget per specimen:

```r
refined <- refine_spatial_labels(xy, labels, samples, workers = 4)
```

The same budget is reused as specimens are processed, so there is no nested
process pool or worker multiplication. Native parallelism is available on macOS,
Linux, and Windows. Changing `workers` does not change labels or inferential
diagnostics; only the `workers` bookkeeping attribute records the effective
budget.

On an 8-core Apple M3, the release benchmark for 100,000 genuinely 3D locations
in four specimens took a median 0.518 seconds with one worker and 0.183 seconds
with four workers (2.83x speedup; three measured repetitions per budget).
Validation, specimen indexing, axis reduction, coordinate normalization, and
native refinement are included; simulation generation is excluded. The script
and raw timings are under `benchmarks/benchmark_release_workers.R` and
`benchmarks/results/release_0.1.0_worker_benchmark/`.

### Genuine 3D coordinates

```r
volume <- simulate_volumetric_domains(
  n = 12000,
  shape = "folded_layers",
  samples = 2,
  seed = 42
)

volume_refined <- refine_spatial_labels(
  volume$xy, volume$labels, volume$samples, workers = 4
)
attr(volume_refined, "dimensions_used")
```

### Independent specimens with overlapping coordinates

```r
xy_one <- sim$xy[sim$samples == levels(sim$samples)[1], , drop = FALSE]
xy <- rbind(xy_one, xy_one)
samples <- rep(c("section_1", "section_2"), each = nrow(xy_one))
labels <- factor(rep(sim$labels[sim$samples == levels(sim$samples)[1]], 2))

joint <- refine_spatial_labels(xy, labels, samples, workers = 4)
attr(joint, "sample_sizes")
```

## Clean a Categorical Mask

For a 2D pixel mask or 3D voxel array, coordinates are generated automatically.
Missing cells remain missing. No source image, intensity, model logit, or clean-mask
training collection is used.

```r
truth <- matrix(rep(c("background", "region"), each = 5000), nrow = 100)
noisy <- corrupt_categorical_mask(
  truth, mechanism = "boundary", rate = 0.15, seed = 4
)
cleaned <- clean_categorical_mask(noisy, workers = 1)

evaluate_mask_cleaning(truth, noisy, cleaned)
```

The evaluator reports multiclass mean IoU, one-grid-step Boundary IoU, rare-class
IoU, correction recall, and damage. Its accuracy obeys the exact identity
`gain = repaired_fraction - damaged_fraction`.

## Diagnostics

The result is an ordinary factor with the input levels. Pointwise diagnostics are
attached without expanding the public interface:

```r
attr(refined, "candidate")
attr(refined, "margin_score")
attr(refined, "required")
attr(refined, "repair_margin")
attr(refined, "atlas_dispersion")
attr(refined, "isolation")
attr(refined, "changed")
```

`margin_score` is an uncalibrated rival-versus-observed evidence contrast, not a
probability or expected loss. `repair_margin >= 0` is the selective change
criterion for a candidate that differs from the observed label.
`atlas_dispersion` is the empirical standard deviation of the deterministic
chart contrasts divided by the square root of the atlas size. It is a
descriptive admission scale, not a calibrated probability or confidence
interval.
`isolation` is the multiclass local-gap protection factor; it is one for the
binary ballot.

The following named summaries are also attached:

```r
attr(refined, "workers")
attr(refined, "dimensions_used")
attr(refined, "labels_changed")
attr(refined, "changed_fraction")
attr(refined, "classes_before")
attr(refined, "classes_after")
attr(refined, "removed_classes")
attr(refined, "sample_sizes")
```

All summaries except `workers` are reported per specimen; class summaries are
named lists. FiberMargin does not force the output to retain every observed
class. `removed_classes` makes any disappearance explicit.

## Evaluation

The package provides planar, gradient, complex-shape, and volumetric simulators,
plus one evaluator used by every benchmark:

```r
evaluate_spatial_refinement(
  truth = sim$truth,
  initial = sim$labels,
  refined = refined,
  boundary = sim$boundary,
  regions = sim$area,
  sparse = sim$sparse,
  method = "FiberMargin"
)
```

`correction_recall` is the fraction of initially wrong labels repaired.
`damage_rate` is the fraction of initially correct labels made wrong. Worst-class,
sparse-region, and boundary accuracy prevent aggregate accuracy from hiding local
failures.

Compact DLPFC and CRC Visium HD coordinate/annotation benchmarks are bundled. CRC
denotes colorectal cancer. DLPFC is redistributed under the `spatialLIBD` terms. The
CRC derivative is CC BY 4.0 and contains only coordinates, 19 author-derived region labels, and stored
corruption recipes; expression counts and tissue imagery are excluded. MERFISH
inputs remain external. See
`inst/extdata/REAL_DATA_LICENSES.md` for sources, attribution, and change notices.

```r
crc <- load_spatial_benchmark(
  "crc", "CRC_random_25_r1", seed = 1040001L
)
mean(crc$labels != crc$truth)
```

## Current Validation

The C++ implementation is evaluated on predefined CRC, DLPFC, and MERFISH
corruption protocols. On the direct-comparator matrices, it obtains 0.8589
accuracy and 0.7806 ARI across 60 CRC corruptions, within 0.0015 and 0.0021 of
multiscale mode. Across 45 DLPFC corruptions it obtains 0.8942 accuracy and
0.8220 ARI. Across 45 MERFISH corruptions it obtains 0.8905 accuracy and
0.8004 ARI, exceeding the strongest fixed comparator by 0.0062 accuracy and
0.0116 ARI. These are controlled label-recovery results, not claims that the
biological annotations are naturally wrong.

The protocol definitions, package-build metadata, per-case ledgers, and matched
plots are under
`benchmarks/results/fibermargin_atlas912_full_roster_v1/`; the source-matched
FiberMargin confirmation is in
`benchmarks/results/fibermargin_atlas912_full_confirmation_v1/`.
Additional simulation and mask panels are retained in the repository and are
re-evaluated only when the production source changes.

See `vignette("fibermargin", package = "fibermargin")` for the workflow,
`vignette("benchmarks", package = "fibermargin")` for benchmark design,
and `vignette("fibermargin_reference", package = "fibermargin")` for the package
reference API.
