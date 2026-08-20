# Human Longevity Bottleneck Model

Minimal reproducibility repository for the npj Aging manuscript:

**Clonal hematopoiesis drives organ bottleneck switching in theoretical human longevity**

## What is included

Only files needed to reproduce or audit the principal analyses reported in the article are included.

- `src/run_efimov_modelII_chip_extension.py`  
  Reconstructs the post-mitotic Efimov Model II brain/heart backbone and runs the frozen CHIP/VAF extension.

- `src/make_figure1.py`  
  Recreates the single five-panel main manuscript figure from the frozen article outputs.

- `config/frozen_config.json`  
  Frozen model parameters and random seed.

- `results/Efimov_literal_reproduction_check.csv`  
  Published-backbone brain, heart and combined reconstruction values.

- `results/ModelIII_MonteCarlo_check.csv`  
  Frozen proliferative-tissue verification values used in panel a and in the Results.

- `results/Efimov_ModelII_CHIP_primary_results.csv`  
  Primary driver-specific T50 and Heart-first/Brain-first probabilities.

- `results/Efimov_ModelII_switching_sweep.csv`  
  Common-HR sweep used to estimate the bottleneck-switching boundary.

- `results/Final_ModelII_primary_and_validation_results.csv`  
  Primary incident-HF, positive CHD and CVD-null triangulation outputs.

- `figures/Figure1_integrated_exact.png`  
  Final main manuscript figure.

## Frozen settings

- Random seed: `20260819`
- Synthetic population: `N = 100000`
- Published non-aging background median: `1759 years`
- Background hazard: `ln(2)/1759`
- CHIP prevalence capped at 30% beyond age 80
- Approximate switching boundary on the exact reconstructed Model II backbone: `HR ≈ 1.27` at `VAF = 10%`

## Reproduce the core Model II plus CHIP analysis

```bash
python src/run_efimov_modelII_chip_extension.py
```

## Recreate the manuscript figure

```bash
python src/make_figure1.py
```

## Interpretation

Theoretical T50 values are model outputs and are not clinical lifespan predictions or estimates of maximum human lifespan.

## Code availability

Repository URL: **https://github.com/soricast/human-longevity-bottleneck-model**

The release used for the manuscript should be tagged `v1.0.0` and may additionally be archived in Zenodo.

## Data availability

No new participant-level clinical data were generated. Empirical inputs were taken from published studies; the repository contains only derived/frozen model outputs needed to reproduce or audit the article.

## License

MIT.
