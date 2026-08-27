# FiberMargin 0.1.0 integration contract

Call `refine_spatial_labels(xy, labels, samples = NULL, workers = NULL)` with one
row per spatial location. `xy` is a finite numeric matrix with exactly two or
three columns, and `labels` contains one non-missing assignment per row.
`samples` may be integer, character, or factor; each observed value defines an
independent coordinate system. Coordinates and class evidence never cross a
sample boundary, even when coordinates overlap numerically.

Within each sample, zero-range coordinate axes are removed. Two varying axes
invoke the 2D operator and three invoke the genuine 3D operator. Fewer than two
varying axes is an error. `workers` is one total cross-platform native CPU budget
for the call. It is reused between samples and does not change results or
diagnostics.

The return value is an ordinary factor in original row order with the input
factor levels and identifiers. Pointwise attributes are `candidate`,
`margin_score`, `required`, `repair_margin`, `atlas_dispersion`, `isolation`, and
`changed`. Summary attributes are `workers`, `dimensions_used`,
`labels_changed`, `changed_fraction`, `classes_before`, `classes_after`,
`removed_classes`, and `sample_sizes`. Except for scalar `workers`, summaries
are named by sample; class summaries are named lists. FiberMargin does not force
class preservation, and any observed class absent after refinement is reported
in `removed_classes`.
