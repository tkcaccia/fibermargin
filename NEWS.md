# fibermargin 0.2.4

- Use `RefineMargin()` for coordinate-and-label refinement. The former
  `refine_spatial_labels()` name remains an exact compatibility alias.
- `RunFiberMargin()` refines an existing clustering in a Seurat, Giotto
  Suite, or SpatialExperiment object and stores the result in a new
  observation metadata column.
- Added an installed vignette with runnable Seurat, SpatialExperiment, and
  Giotto Suite examples.
- Added a public API benchmark covering simulated 2D and 3D specimens and
  the bundled DLPFC and CRC data. No changes to native refinement or scoring.

# fibermargin 0.2.3

* Evaluate masks on the fixed reference/input footprint, rejecting missing
  cleaned labels on that footprint. Reject observed missing factor levels
  in spatial benchmarks and record invalid predictions as method failures.
* Seeded simulations, corruptions, and benchmark runs restore the caller's
  random-number state without changing the refinement algorithm.

# fibermargin 0.2.2

* Generate the same Hilbert keys with specialized 2D and 3D code, and select
  the exact normalization quantiles without sorting entire axes.
* Reuse one transport buffer per worker, build each path on demand, and compute
  dispersion only for the selected class. The public API and the refinement
  rule are unchanged.

# fibermargin 0.2.1

* Reject observed missing factor levels in labels and specimen identifiers.
  Preserve unused input label levels in the returned factor.
* Normalize without an absolute tolerance in the input coordinate units,
  including for small spans and coordinates close to the double-precision
  range limit. Use the full axis range if both robust quantiles coincide.

# fibermargin 0.2.0

* Resolve Hilbert key ties by coordinates. Coincident multiclass sites share
  peer evidence symmetrically without counting the query itself.
* Give binary neighbourhoods a deterministic ordering at equal distances.
* Index rows by specimen once, and join started CPU workers even if a later
  thread launch fails. The maximum supported worker budget is 64.
* Avoid starting additional threads for short admission passes.

# fibermargin 0.1.0

* Stabilized the four-argument `refine_spatial_labels()` public API.
* Added strict validation for coordinate dimensionality, finite values, labels,
  specimen identifiers, worker budgets, and per-specimen coordinate variation.
* Made coordinate dimensionality specimen-specific. Constant axes are removed
  independently, while genuinely three-dimensional specimens retain all three
  coordinates.
* Made specimen isolation explicit in the native engine, including local class
  encodings and support for completely overlapping coordinate systems.
* Replaced platform-specific specimen forking with one cross-platform C++ thread
  budget. Results and diagnostics are deterministic across worker counts.
* Added pointwise and per-specimen summary diagnostics, including class-loss
  reporting without imposing class preservation.
* Added 2D, genuine 3D, multi-specimen, invariance, memory, and parallelism tests
  and expanded the README, vignette, reference documentation, and manuscript.
* Removed persistent vignette option changes and replaced the internal fixed
  CRC benchmark seed with the optional user-controlled `seed` argument.
