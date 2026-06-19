# A covariate's effect on a focal taxon can have opposite signs in a Beta
# (relative-abundance) model versus a CLR model, for the same data. The two
# models measure change against different references: the Beta model against
# the taxon's share of the arithmetic total, the CLR model against the community
# geometric mean. When rare and abundant taxa respond differently to a
# covariate, those references move in opposite directions and the focal
# coefficient flips sign, even while the taxon's absolute abundance is rising.
#
# Reproduces the example in Appendix S3 (embedded in manuscript/manuscript.tex).

library(betareg)
set.seed(1)

n <- 4000
x <- rnorm(n)                          # standardized covariate (temperature)

# True process: each taxon's log abundance is linear in x. The focal taxon
# increases moderately, a rare taxon increases strongly, an abundant taxon flat.
b  <- c(focal = 1, rare = 4, abundant = 0)        # response slopes
c0 <- c(focal = 0, rare = -4, abundant = 0)       # baseline log-abundance

a <- exp(sapply(1:3, function(k) c0[k] + b[k] * x + rnorm(n, 0, 0.05)))
p <- a / rowSums(a)                    # relative abundances (sum to one)

# Fit the focal taxon under each model.
beta_coef <- coef(betareg(p[, 1] ~ x | 1, link = "cloglog"))[["x"]]
clr       <- log(p) - rowMeans(log(p))
clr_coef  <- coef(lm(clr[, 1] ~ x))[["x"]]

# Exact slopes for this log-linear process: focal slope minus a reference.
p_bar <- colMeans(p)
results <- data.frame(
  model    = c("Beta (cloglog)", "CLR"),
  fitted   = c(beta_coef, clr_coef),
  exact    = c(b[1] - sum(p_bar * b),   # abundance-weighted reference
               b[1] - mean(b)),         # geometric-mean reference
  reference = c("arithmetic total", "geometric mean")
)
print(results, row.names = FALSE, digits = 3)
