# fibermargin 0.2

- `FiberMargin()` refines spatial labels from 2D or 3D coordinates and
  reports candidate labels, edit margins, and specimen-level diagnostics.
- `RunFiberMargin()` writes refined labels to Seurat, Giotto, and
  SpatialExperiment objects. `clean_categorical_mask()` supports regular
  masks and volumes.
- Multiclass refinement offers standard and short transport presets.
  The binary ballot is the same for both presets.
- Specimens are refined independently, with constant axes removed within
  each specimen. Labels and pointwise diagnostics are deterministic across
  CPU worker counts.
- The package includes simulation, scoring, and benchmarking tools and
  coordinate-and-label datasets for DLPFC and colorectal cancer.
- Vignettes cover coordinate, spatial-object, and benchmark workflows.
