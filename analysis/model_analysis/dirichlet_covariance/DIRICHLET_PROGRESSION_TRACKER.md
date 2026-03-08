# Dirichlet Model Progression Tracker

## Current Status: ✅ WORKING BASE MODEL
- **File**: `01_fitModels_dirichlet_real_data.r`
- **Status**: MCMC completes successfully with real data
- **Dimensions**: 10 plots × 9 dates × 4 species (121 observations)
- **Samplers**: 192 (manageable)
- **Samples**: 100 collected successfully

## Progression Plan

### Phase 1: ✅ COMPLETED - Basic Working Model with Convergence
- [x] Real data loading and basic filtering
- [x] Simple timepoint/plot/site mapping
- [x] Top 3 taxa + "other" category
- [x] Basic Dirichlet model structure
- [x] Simple priors
- [x] MCMC completion
- [x] **CONVERGENCE TEST**: 3 chains, 1000 iterations each
  - **Mean GBR**: 1.003 (excellent, well below 2.0 threshold)
  - **All individual GBR**: < 1.02 (excellent convergence)
  - **Mean ESS**: 1,476 (very good effective sample size)
  - **Min ESS**: 1,223 (good minimum effective sample size)

### Phase 2: 🔄 IN PROGRESS - Add Environmental Covariates
- [ ] Replace synthetic environmental data with real predictor data
- [ ] Add proper environmental covariate loading
- [ ] Test with real environmental predictors
- [ ] Verify MCMC still completes

### Phase 3: 🔄 PENDING - Add Complex Model Structure
- [ ] Add proper alpha dynamics (temporal evolution)
- [ ] Add environmental effects on alpha
- [ ] Add site effects
- [ ] Test with more complex model

### Phase 4: 🔄 PENDING - Add Convergence Features
- [ ] Add convergence-based sampling
- [ ] Add proper error handling
- [ ] Add parallel execution capability
- [ ] Test with full framework

### Phase 5: 🔄 PENDING - Scale Up Data
- [ ] Increase number of plots (10 → 50)
- [ ] Increase number of dates (9 → 20)
- [ ] Test with larger dataset
- [ ] Optimize for performance

### Phase 6: 🔄 PENDING - Full Integration
- [ ] Integrate with main script framework
- [ ] Add proper output saving
- [ ] Add logging and monitoring
- [ ] Final testing

## Testing Protocol
Each phase will be tested with:
1. MCMC completion (no hanging)
2. Sample collection verification
3. Parameter monitoring confirmation
4. Performance metrics (sampler count, runtime)

## Current Model Structure
```r
# Simple Dirichlet model
for (i in 1:N.core) {
  y[i, 1:N.spp] ~ ddirch(alpha[plot_num[i], 1:N.spp, timepoint[i]])
}

# Simple alpha generation
for (s in 1:N.spp) {
  for (p in 1:N.plot) {
    for (t in 1:N.date) {
      alpha[p, s, t] ~ dgamma(1, 1)
    }
  }
}

# Simple priors
for (s in 1:N.spp) {
  sigma[s] ~ dgamma(1, 1)
  intercept[s] ~ dnorm(0, 1)
  beta[s, 1:4] ~ dmnorm(zeros[1:4], omega[1:4, 1:4])
}
```

## Target Model Structure (from 01_fitModels_betaReg.R)
```r
# Complex Dirichlet model with environmental effects
for (i in 1:N.core) {
  y[i, 1:N.spp] ~ ddirch(alpha[plot_num[i], 1:N.spp, timepoint[i]])
}

# Plot-level process model with environmental effects
for (s in 1:N.spp) {
  for (p in 1:N.plot) {
    # Initial condition
    alpha[p, s, plot_start[p]] ~ dgamma(1, 1)
    
    # Dynamic evolution with environmental effects
    for (t in (plot_start[p] + 1):N.date) {
      log(Ex[p, s, t]) <- rho[s] * log(max(alpha[p, s, t - 1], 0.001)) +
        beta[s, 1] * temp[plot_site_num[p], t] +
        beta[s, 2] * mois[plot_site_num[p], t] +
        beta[s, 3] * pH[p, t] +
        beta[s, 4] * pC[p, t] +
        beta[s, 5] * relEM[p, t] +
        beta[s, 6] * LAI[plot_site_num[p], t] +
        site_effect[plot_site_num[p], s] +
        intercept[s]
      
      alpha[p, s, t] ~ T(dnorm(mean = Ex[p, s, t], sigma[s]), 0.001, Inf)
    }
  }
}

# Complex priors
for (s in 1:N.spp) {
  sigma[s] ~ dgamma(2, 15)
  intercept[s] ~ dt(0, 0.3, df = 3)
  rho[s] ~ dbeta(2, 2)
  beta[s, 1:8] ~ dmvt(zeros[1:8], omega[1:8, 1:8], df = 3)
}

# Site effects
sig ~ dgamma(2, 8)
for (s in 1:N.spp) {
  for (k in 1:N.site) {
    site_effect[k, s] ~ dt(0, sig, df = 3)
  }
}
```

## Notes
- Each phase builds on the previous one
- Testing is done after each phase
- If a phase fails, we revert and try a simpler approach
- Goal is to reach full complexity while maintaining MCMC completion
