## Update

VsusP 1.1.0 updates version 1.0.0. The changes include corrected S2M stopping
counts, additional selection summaries, Gaussian fitting and diagnostic plots.
The maintainer and license are unchanged.

## Test environment

* macOS 14.6 arm64, R 4.5.0: `R CMD check --as-cran` on the final source archive,
  with online incoming and system-clock checks enabled.

## Check results

0 errors, 0 warnings, 1 note.

```text
checking for future file timestamps ... NOTE
unable to verify current time
```

R could not verify the system clock through its external time services. No
future-dated package files were reported. CRAN incoming feasibility, installed
tests (578 assertions; no failures, warnings or skips), PDF/HTML documentation
and rebuilt vignette checks passed.

## Reverse dependencies

The CRAN source index lists no packages with VsusP in Depends, Imports,
LinkingTo, Suggests or Enhances. No reverse-dependency package checks apply.
