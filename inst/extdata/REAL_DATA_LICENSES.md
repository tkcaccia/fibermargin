# Real benchmark data and redistribution

## Bundled DLPFC benchmark

`dlpfc_benchmark` contains only coordinates, manual cortical-layer labels,
sample identifiers, derived boundary flags, and FiberMargin corruption scenarios.
It contains no expression matrix or histology image. The source is the
`spatialLIBD` Bioconductor data package, licensed under Artistic-2.0:

- https://bioconductor.org/packages/spatialLIBD
- Maynard et al. (2021), DOI: 10.1038/s41593-020-00787-0
- Pardo et al. (2022), DOI: 10.1186/s12864-022-08601-w

The complete upstream license is included as `LICENSE-spatialLIBD` in this
directory.

The FiberMargin corruption matrices and derived boundary/sparsity flags are
distributed under the package license.

## Bundled colorectal benchmark

`colorectal_benchmark` is a compact derivative of:

> Visium HD Spatial Gene Expression Library, Human Colorectal Cancer (FFPE),
> HD Spatial Gene Expression dataset analyzed using Space Ranger v4.0.1,
> 10x Genomics (2025, July 3).

- Source: https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-human-crc-v4
- Citation guidance: https://www.10xgenomics.com/support/software/cell-ranger/latest/miscellaneous/cr-citations
- License: Creative Commons Attribution 4.0 International (CC BY 4.0),
  https://creativecommons.org/licenses/by/4.0/

The derivative retains coordinates and 19 author-derived WSI region labels for
194,541 annotated locations. It excludes expression counts, unannotated locations,
and tissue imagery, and adds boundary, adjacency, sparse-class, and deterministic
corruption metadata. The author-derived labels and compact derivative are released
under CC BY 4.0. See `ATTRIBUTION-colorectal.md` for the complete attribution and
change notice. No endorsement by 10x Genomics is implied.

## Data not bundled

The Moffitt et al. MERFISH measurements deposited at Dryad are CC0
(DOI: 10.5061/dryad.8t8s248). The publication benchmark, however, uses a
processed anatomical-domain annotation obtained from the BASS analysis
repository, which does not state an explicit redistribution license. Those
processed labels are therefore not copied into the package.
