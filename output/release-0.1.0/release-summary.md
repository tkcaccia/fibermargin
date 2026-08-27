# FiberMargin 0.1.0 release summary

Release date: 2026-08-27

## Package verification

- Package version: 0.1.0
- Unit and integration tests: 32 test blocks, 243 expectations, 243 passed
- `R CMD check --as-cran`: 0 errors, 0 warnings, 1 environment note
- The note reports `ResultBundle_*.xcresult` and `xcrun_db` files left in the
  temporary directory by the Apple compiler toolchain. Package installation,
  examples, tests, vignettes, manual generation, and compiled-code checks pass.

## Worker benchmark

The benchmark uses 100,000 genuinely three-dimensional points distributed over
four independent specimens. It includes input validation, specimen indexing,
per-specimen axis reduction and normalization, native refinement, and diagnostic
assembly. Simulation time is excluded. Three measured repetitions follow one
warm-up execution.

| Total worker budget | Median elapsed time (seconds) |
|:--|--:|
| 1 | 0.518 |
| 4 | 0.183 |

Median speedup: 2.83x.

## Public contract

`refine_spatial_labels(xy, labels, samples = NULL, workers = NULL)` accepts a
finite numeric coordinate matrix with two or three columns. `samples` partitions
rows into independent coordinate systems. Constant coordinate axes are removed
separately within each sample, and every sample must retain at least two varying
axes. Genuine 3D specimens use all three coordinates. `workers` is one total,
deterministic CPU budget.

The return value is an ordinary factor in original row order with the input factor
levels and row identifiers. Pointwise attributes are `candidate`, `margin_score`,
`required`, `repair_margin`, `atlas_dispersion`, `isolation`, and `changed`.
Summary attributes report the worker budget, dimensions used, changes, class sets,
removed classes, and sample sizes. The method does not force every input class to
survive.
