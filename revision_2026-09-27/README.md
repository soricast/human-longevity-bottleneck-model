# Reviewer revision computational package

This directory contains only code and derived outputs for the 2026-09-27 reviewer revision. Manuscript, response letter, and journal Supplementary Data are submitted separately to the journal. The original v1.0.0 release remains historically intact.

## Files

- `src/`: original primary Model II extension and all new simulation and figure scripts.
- `results/`: generated numerical CSV outputs and cached intrinsic curves used for the revised analysis.
- `figures/`: main Figures 1–5 and Supplementary Figure S1 (`Figure6_VAF_population.png`) and S2 (`Figure7_KV_extrapolation.png`).
- `METHODS_AND_PROVENANCE.md`: detailed assumptions, provenance, script order, and interpretation limits.
- `SHA256SUMS.txt`: checksums of the files in this directory.

## Reproduce

Run from this directory in an environment installed with `pip install -r requirements.txt`:

```bash
python src/run_replicates.py
python src/make_new_figures.py
python src/run_composition_multiclone.py
python src/run_joint_uncertainty.py
python src/run_correlation_sensitivity.py
python src/make_uncertainty_figure.py
python src/make_vaf_k_figures.py
```

Scripts write regenerated outputs to `analysis_outputs/`; frozen revision outputs are in `results/`. The ten-seed ranges quantify Monte Carlo variation at fixed model inputs; joint-parameter quantiles depend on published marginal intervals and modeling assumptions. Non-cardiac coefficients are illustrative, and the non-DNMT3A cardiac HR is a composite proxy. No participant-level data are included.
