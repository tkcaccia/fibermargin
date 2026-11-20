## Update to fibermargin 0.3.0

The default algorithm and outputs are unchanged. The three refinement
interfaces now accept an optional short transport preset for multiclass
coordinates. The package retains one implementation; the preset changes
only a fixed decay range. No new mandatory dependencies or data downloads
were added. The vignettes and reference pages document the choice.

## R CMD check

Local source tarball checked with `R CMD check --as-cran` on macOS
arm64 under R 4.6.0: 0 errors, 0 warnings, 1 note.

The incoming feasibility note listed two items:

* GiottoClass is an optional enhancement distributed by Giotto Suite.
  Core functionality, installation and vignette construction do not
  require it. Integration examples and tests detect its availability.
* The 10x Genomics dataset URL in `man/colorectal_benchmark.Rd` returned
  HTTP 429 to the automated URL check. It is the original public data
  attribution URL. The included coordinate-and-label derivative has
  source, license, and attribution information in `inst/extdata`.
