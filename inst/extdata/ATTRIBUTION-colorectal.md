# Colorectal Visium HD benchmark attribution

## Source

Visium HD Spatial Gene Expression Library, Human Colorectal Cancer (FFPE),
HD Spatial Gene Expression dataset analyzed using Space Ranger v4.0.1,
10x Genomics (2025, July 3).

Source: <https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-human-crc-v4>

10x Genomics citation guidance:
<https://www.10xgenomics.com/support/software/cell-ranger/latest/miscellaneous/cr-citations>

## License

The source dataset is licensed under the Creative Commons Attribution 4.0
International license (CC BY 4.0):
<https://creativecommons.org/licenses/by/4.0/>.

The 19 WSI annotation labels added by the project authors and the compact package
derivative are also distributed under CC BY 4.0. This attribution does not imply
endorsement by 10x Genomics.

## Changes

The package derivative retains spatial coordinates and 19 author-derived WSI region
labels for annotated locations. Unannotated locations, expression counts, and tissue
imagery are excluded. Boundary, adjacent-region, sparse-class, and deterministic
corruption metadata were computed for label-refinement evaluation. The benchmark
therefore must not be treated as an unmodified 10x Genomics data release.
