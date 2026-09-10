## Resubmission

This update addresses the requested CRAN changes:

* No reference was added to `DESCRIPTION` because no publication currently
  describes the FiberMargin method.
* Vignettes no longer change `options()`, and every graphical `par()` change is
  paired with restoration of the complete prior graphical state.
* `load_spatial_benchmark()` no longer sets an internal fixed seed. Its optional
  `seed = NULL` argument leaves seed selection to the user; a supplied integer
  provides reproducible CRC corruption generation.

## R CMD check results

Tested with R 4.6.0 on macOS arm64 using:

```text
R CMD check --as-cran fibermargin_0.1.0.tar.gz
```

Result: 0 errors, 0 warnings, 1 note.

The incoming URL check received HTTP 429 (Too Many Requests) from the 10x
Genomics dataset page documented in `man/colorectal_benchmark.Rd`:

```text
https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-human-crc-v4
```

The URL is the public source page for the bundled CRC benchmark derivative.
The remote server rate-limited the automated check; all package checks
completed successfully.

## Package scope

FiberMargin contains a C++17 implementation and an R interface. The source
tarball includes two compact coordinate-and-label benchmark derivatives. Their
source, attribution, and redistribution terms are recorded in `inst/extdata`.
No external data are downloaded during installation, examples, tests, or
vignette construction.
