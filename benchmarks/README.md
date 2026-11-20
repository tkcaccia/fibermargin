# FiberMargin benchmarks

Install the package from the source revision being audited before running
these commands from the repository root:

```sh
Rscript benchmarks/benchmark_public_api.R /tmp/fibermargin-api-benchmark
Rscript benchmarks/benchmark_release_workers.R /tmp/fibermargin-worker-benchmark
```

The public API audit runs simulations in 2D and genuine 3D, including a
thin-sheet fixture, and one packaged corruption each for DLPFC and CRC.
It scores the standard and short presets against the same input. It checks
the compatibility alias for the default and verifies exact equality between
one and four workers for both presets, apart from the worker-count attribute.
The report uses the default public benchmark scorer and writes all metrics
and R session information. These are regression checks, not new estimates of
real-world annotation accuracy.

The worker benchmark uses the same 100,000-location simulated volume as the
archived release timing, with a warm-up and three measured calls per worker
budget. To compare binary, multiclass, and multi-specimen predictions with an
older installed release, use the separate-library instructions in
`benchmark_native_performance.R`.

Earlier experiments under `marginsvm_*` and their result ledgers are historical
and may require the associated archived software. The scripts above exercise
the current package API; they do not overwrite archived results.
