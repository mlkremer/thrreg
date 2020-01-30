# thrreg 0.1.0 (January 29, 2020) 

Updates and changes in functions compared to 
[original source code](https://www.ssc.wisc.edu/~bhansen/progs/ecnmt_00.html) 
by Bruce E. Hansen.

### `thr_test_hom()` and `thr_test_het()`
* New arguments: `var.names`, `cr`, `graph`, `quick`
* A `data.frame` can be passed
* Arguments `yi`, `xi`, `qi` may be of type `integer` or `character`
* Output in LaTeX format 

### `thr_est()`
* New arguments: `test.pvalue`, `conf2`, `nonpar`, `graph`, `signif.level`, 
  `digits`, `integer.digits`, `header`, `output.short`, `signif.legend`
* A `data.frame` can be passed
* Arguments `yi`, `xi`, `qi` may be of type `integer` or `character`
* Loop over different significance levels `c(.9, .95, .99)` for estimation of 
  coefficients
* Definition of notation for significance levels `c(.9, .95, .99)` for 
  threshold test and estimated coefficients
  * for `"stars"`: number of stars
  * for `"colors"`: red/blue tones for positive/negative estimates
* Constructs a summarizing threshold regression table with estimates, 
  standard errors, and significance levels for both regimes
* Changed `"Constant"` to `"Const"` for shorter output in tables
* Output in LaTeX format
