# RSQ vs CRPS discrepancy for modeled vs random site effects

## Observation

- **RSQ** (fig_newtime_newsite.r → "RSQ by Site Prediction Type"): Modeled site effects show higher RSQ than random site effects.
- **CRPS** (fig_site_effect_improvement.r): CRPS improvement from modeled vs random site effects is often zero or negative (no improvement or worse).

## Data source

Both metrics are computed on the **same hindcasts** in `08_calculateScoringMetrics.r`: for each `(model_id, site_prediction)` the same `(truth, mean, med, sd)` are passed to `robust_add_scoring_metrics()`, which returns both RSQ and CRPS (and others). So the discrepancy is not from different data or aggregation units.

## Why the metrics can diverge

1. **What each metric uses**
   - **RSQ** (in `robust_add_scoring_metrics.R`): `RSQ = summary(lm(observed ~ mean_predicted))$r.squared`. With `use_median = TRUE`, the point forecast used for RSQ (and RMSE, etc.) is the **median** of the predictive distribution. RSQ measures how much variance in observations is explained by this **point forecast** only.
   - **CRPS**: Uses the **full predictive distribution** (mean and sd) via `crps_norm(observed, mean_predicted, sd_predicted)`. CRPS rewards accurate **and** well-calibrated probabilistic forecasts (location and spread).

2. **Interpretation**
   - Modeled site effects can yield **better point forecasts** (medians closer to truth on average) → **higher RSQ** than random.
   - The **full predictive distribution** for modeled site effects can still be **poorly calibrated**: e.g. too narrow (overconfident) or biased spread. CRPS penalizes bad spread and location. Random site effects often use a simpler, wider spread (site mean ± uncertainty), which can score better on CRPS when the modeled forecast is overconfident.
   - So: better point accuracy (RSQ) does **not** imply better probabilistic performance (CRPS) if uncertainty is mis-specified.

3. **Summary**
   - RSQ reflects only point-forecast skill; CRPS reflects full distributional skill.
   - The figures are consistent: modeled site effects add some point-forecast value (RSQ) but do not improve (or can worsen) probabilistic skill (CRPS), likely due to calibration/spread of the predictive distribution.

## References in code

- RSQ/CRPS definition: `analysis/model_analysis/robust_add_scoring_metrics.R` (RSQ line 79, CRPS lines 54–67).
- Scoring by `site_prediction`: `analysis/model_analysis/08_calculateScoringMetrics.r` (grouping by `site_prediction` ~line 216).
- RSQ figure: `analysis/create_figs/fig_newtime_newsite.r` (rank_vals from `scoring_metrics`, RSQ by `site_prediction`).
- CRPS improvement figure: `analysis/create_figs/fig_site_effect_improvement.r` (CRPS from `scoring_metrics_cv`, improvement = (Random - Modeled)/Random).
