# Package checks

Run `check.R` against a built source archive and an empty results directory:

```sh
Rscript tools/check.R VsusP_1.1.0.tar.gz quick results
```

Use the actual archive filename. `quick` checks the installed package and tests
without TeX or online incoming checks. `docs` also rebuilds the vignette and
checks PDF/HTML documentation. `full` adds online incoming and system-clock
checks. Results retain logs, package versions, archive identity and test counts.

CI uses Windows, macOS and three Linux R versions. Linux R release also
checks documentation; other jobs avoid TeX. `check-results.R` rejects incomplete
checks, warnings, missing tests, skipped tests and zero executed assertions.
Coverage runs only for relevant source/test changes and retains its XML artifact;
Codecov upload is separate and requires the configured token.
