# Paired simulation extension

Run from this directory with Python and NumPy, SciPy, pandas, matplotlib and openpyxl installed:

```bash
python src/run_replicates.py
python src/make_new_figures.py
```

The original Model II primary script is included in `src/`. Ten independent synthetic populations use seeds 20260819–20260828, each with 100,000 individuals. Within each seed, all scenarios reuse the same carrier status, driver, onset, growth and three exponential failure thresholds. The numerical Model II organ survival curves are held fixed across seeds and generated with the original fixed organ-curve seed. The time grid is 0.25 years through 1,000 years, matching the original script. The first replicate reproduces its primary T50 and pairwise heart-first probabilities.

`null` removes CHIP cardiac coupling. `cardiac` applies the original clone-dependent cardiac multiplier. `brain_half`, `brain_equal` and `brain_double` also apply the multiplier to brain hazard with exponents 0.5, 1 and 2. `background_equal` adds an unrelated cumulative hazard equal to the integrated cardiac excess hazard, while retaining the cardiac multiplier. These non-cardiac coefficients are hypothetical stress tests and are not estimated from participant data.

`P_heart_before_brain` compares theoretical organ times even when background failure comes earlier. `P_terminal_heart` counts the heart as the earliest among heart, brain and background. `delta_*_matched` subtracts the null result for the same seed and assigned driver group. The 10-seed minimum and maximum describe simulation variation at fixed model inputs; they are not confidence intervals for biological parameters. Baseline organ-curve Monte Carlo uncertainty is not varied.

Age-specific `ratio_heart_brain` is the discrete Model II heart hazard increment divided by the brain increment, multiplied by the mean cardiac multiplier among all assigned carrier trajectories at each age. It is not survivor-weighted or an observed clinical cause-specific hazard ratio. The original single-seed Weibull intrinsic-hazard sensitivity data are used only for Figure 4 and remain a separate model representation.

Results in `results/` and figures in `figures/` correspond to this revision. The frozen published analysis and separate sensitivity workbook are retained in the submitted supplementary files.

The composition and multiple-clone extension is run with `python src/run_composition_multiclone.py`. ARIC follow-up clone occurrence counts (540/258/97/18) come from Uddin et al. 2024 Supplementary Data 2; 383/457 baseline carriers had one clone. Counts are renormalized among these four genes only. The alternative multiple-clone rules use the maximum or product of the two clone cardiac multipliers and are stress tests. This extension draws second-clone variables before the failure thresholds, so its nominal single-clone Lothian realization is a new matched simulation and need not equal the original 20260819 frozen run. The original publication source is https://doi.org/10.1038/s41467-024-52302-9 .

Joint uncertainty: run `python src/run_joint_uncertainty.py` then `python src/run_correlation_sensitivity.py` and `python src/make_uncertainty_figure.py`. This uses 200 independent lognormal parameter combinations from the reported marginal brain/heart mutation-rate and adjusted heart-failure HR intervals, a shared 50,000-person synthetic population, and a 0.5-year grid. Organs' cumulative hazards are time-rescaled by sampled mutation-rate ratios; this corresponds to a uniform scaling of their original lognormal loss-rate distributions. The association multiplier is compared when continued indefinitely and when set to one after age 100. The latter is a bounding assumption. Correlation variants use illustrative brain/heart log-rate correlations of -0.5 and +0.5, with 60 draws each. These are conditional model quantiles, not causal or population confidence intervals. Parameter draws are separate from the ten-seed fixed-parameter analysis.

Reviewer 3 visualization: run `python src/make_vaf_k_figures.py` to generate driver-specific VAF median/10th–90th percentile bands and the anchored K_V cardiac multiplier curves. The representative trajectory and distribution are based on synthetic carrier assignments, not observed VAF distributions. K_V=1 percentage point is an illustrative curvature choice; all plotted curves reproduce HR=1.52 at 10% VAF. The original one-at-a-time Weibull K_V sensitivity rows are included separately in `results/KV_surrogate_sensitivity.csv`.
