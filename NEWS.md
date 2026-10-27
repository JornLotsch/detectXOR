# detectXOR 0.2.0

## Changes to the detection rule

* The chi-square comparison of observed and expected tile counts and the fixed
  absolute-difference threshold (`abs_diff_threshold`) were replaced. A pair is
  now tested with the XOR purity: the balanced accuracy of the XOR rule within
  the four critical tiles, maximized over both label assignments and both XOR
  manifestations (default, rotated). A pair is called XOR-positive if
  - the purity exceeds chance in a class-label permutation test, adjusted for
    the number of variable pairs (Westfall-Young max-T by default,
    `p_adjust_method`, `n_perm`);
  - the one-sided lower confidence bound of the purity reaches
    `purity_threshold` (default 0.8);
  - each class populates both of its expected tiles (`min_tile_share`, default 0.05);
  - the class-wise Kendall tau criterion is met (`check_tau`, as before);
  - optionally, all four group-wise Wilcoxon tests are significant
    (`compute_axes_parallel_significance`, now FALSE by default). Previously
    these tests were computed but did not affect the decision.
* `abs_diff_threshold` is deprecated and ignored.
* The results table has new columns (`xor_purity`, `purity_lower`,
  `coverage_ok`, `p_value`, `p_adjusted`, `candidate_pattern`,
  `axis_tests_passed`); `chi_sq` and `chi_sq_p_value` were removed.

## Bug fixes

* Detection no longer depends on which class label is coded first.
* Class columns not named `"class"` no longer cause an error.
* `extreme_handling = "remove"` is accepted as a synonym of `"trim"`, and
  trimming no longer considers the class column.
* Tile assignment for the Wilcoxon tests now uses the same grid as the main
  analysis; pairs of binary variables are exempt from these tests.
* The sequential fallback of the parallel processing wrapper was repaired.
* Winsorization handles missing values when `use_complete = FALSE`.

## Other

* New argument `split_quantiles` to change the tile boundaries (default: tertiles).
* New argument `balance_classes` to down-sample the larger class before screening,
  recommended for marked class imbalance.
* Unit tests were added (testthat).
* Examples use `XOR_data` instead of the three-class `iris` data.
