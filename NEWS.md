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
