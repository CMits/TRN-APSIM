# R side of the reproduction

Verified on **R 4.6.1** (Windows 11).

## What is actually required

**Base R only.** `scripts/run_psoup_sweep.R` sources
`02_networks/Hand_Curated/psoup_r/model/psoup_driver.R` and compiles the committed
equations in `model/model_baseline/nextStep.R`. Neither file loads a package, so a
plain R install runs the hand-curated arm.

## Optional: the PSoup package

Only needed to **regenerate** `model_baseline/nextStep.R` from the SBGN diagram
(`02_networks/Hand_Curated/Dataset S18_Benchmark_Model_Original.sbgn`) via
`model/build_network.R`, which calls `PSoup::convertSBGNdiagram()`. This pipeline
does not do that — it runs the committed equations — so the package is not
required to reproduce any result here.

```r
install.packages(c("igraph", "xml2"))
remotes::install_github("NicoleZFortuna/PSoup",
                        ref = "e74c2b79dc19c3da63e38311b1bc1821949d55ae")
```

`run_all.ps1` phase 00 reports whether the package is present but does not fail
without it.
