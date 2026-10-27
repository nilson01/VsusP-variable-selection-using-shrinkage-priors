# Selection pilot

This exploratory check uses 60 generated Gaussian datasets: 10 for each
combination of two signal coefficients `(0,0)`, `(0.25,-0.25)` or `(1,-1)` and
AR(1) predictor correlation `rho = 0` or `0.8`. Each dataset has 80 observations
and 10 predictors, standardized using sample standard deviations before fitting.
The other eight coefficients are zero; residual standard deviation is one.

Horseshoe fits use `bayesreg` 1.3, one worker, 1,000 retained draws, burn-in 500
and thinning 5. SAVS uses posterior means. S2M comparisons use the same draws
at the fixed illustrative thresholds `0.1`, `0.5` and `1`. These thresholds
were not optimized. Results from one dataset are paired, not independent runs.

| File | Contents |
|---|---|
| `conditions.csv` | Six generating conditions |
| `pilot-status.csv` | All attempted cases, seeds, warnings, failures and timing |
| `pilot-results.csv` | Dataset-level outcomes for every selector |
| `summary.csv` | Means, Monte Carlo standard errors and 95% Wilson intervals |
| `run-selection-pilot.R` | Reproduce all 60 cases in a new output directory |
| `summarize-selection-pilot.R` | Aggregate complete results |
| `plot-selection-pilot.R` | Draw the figure from the summary table |

From a source checkout with VsusP installed:

```sh
Rscript inst/validation/run-selection-pilot.R /path/to/new-output
```

The runner refuses to reuse an occupied directory. It saves data, posterior
draws and every failed case. Summary generation stops if any case fails; no
failed case is replaced or silently excluded. The recorded run used R 4.5.0
on macOS arm64. Small numerical differences may arise on other platforms.

## Interpretation

SAVS selected no variables in all 20 null datasets. Weak signals were often
missed, including both coefficients in every correlated weak-signal case.
For correlated strong signals, SAVS recovered the exact support in 8 of 10
cases. Results from 10 datasets per condition are imprecise; 10/10 successes
has a 95% Wilson interval of approximately 0.72 to 1.00.

The outcome table records false positives, false negatives, false discovery
proportions, empty selections and exact-support recovery. False discovery
proportion is zero for an empty selection. Monte Carlo standard errors use
independent datasets, not individual predictors. Posterior-mean selection
from the first and second halves disagreed in 2 of 60 cases; agreement is not
a convergence certificate. This short single-chain pilot does not establish
false-discovery control, posterior convergence or general performance.

For the unit check, alternate predictor columns were multiplied by `0.1` and
`10`, and their estimates divided by the same factors. The fitted linear
predictor is unchanged, but raw SAVS selections can change. Undoing the
transformation restored all original selections. This is a post-processing
check, not evidence that fitting with different prior scales is invariant.

SAVS follows Ray and Bhattacharya (2018), equation (2) and Algorithm 1:
<https://arxiv.org/abs/1810.09004>. S2M follows Li and Pati (2017):
<https://doi.org/10.1016/j.csda.2016.10.008>. Computational options and their
limits are described in the function help.
