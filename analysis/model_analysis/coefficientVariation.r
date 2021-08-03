# Looking at coefficients of variation (volatility) for microbial time series

install.packages("cvequality")
library(cvequality)

# Read in microbial abundances
d <- c(readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_16S_2021.rds"), 
			 readRDS("/projectnb/talbot-lab-data/zrwerbin/temporal_forecast/data/clean/cal_groupAbundances_ITS_2021.rds"))

# Only keep functional groups; subset to one rank.
ranks.keep <- names(d)
ranks.keep <- ranks.keep[!grepl("bac|fun$", ranks.keep)]
rank.name <- ranks.keep[1]
rank.df <- d[[rank.name]] 
cellulolytic <- d$cellulolytic$cellulolytic
acidobacteriota <- d$phylum_bac$acidobacteriota

# why aren't these working??
# if you are reading this. yes you. tell me what I am doing wrong here?
asymptotic_test(seed = 1, cellulolytic, acidobacteriota)
mslr_test(nr = 10000, x = cellulolytic, y = acidobacteriota)

# tibble::glimpse(cellulolytic)
# num [1:2810] 0.0521 0.0587 0.0534 0.0392 0.0853 ... # all numeric, no NAs
# 
# tibble::glimpse(acidobacteriota)
# num [1:2810] 0.189 0.179 0.162 0.163 0.175 ... # all numeric, no NAs

# cvequality::mslr_test() function code

x <- cellulolytic
y <- acidobacteriota

if (!is.numeric(x) && !is.numeric(y) && !is.character(y)) {
	warning("x is not numeric or y is not numeric or character: returning NA")
	return(NA_real_)
}
if (anyNA(x)) {
	warning("x cannot contain any NA values: returning NA")
	return(NA_real_)
}
if (anyNA(y)) {
	warning("y cannot contain any NA values: returning NA")
	return(NA_real_)
}
n <- data.frame(table(y))$Freq
s <- aggregate(x, by = list(y), FUN = sd)$x
x <- aggregate(x, by = list(y), FUN = mean)$x
k <- length(x)
gv <- as.vector(nr)
df <- n - 1
xst0 <- LRT_STAT(n, x, s)
uh0 <- xst0[1:k]
tauh0 <- xst0[k + 1]
stat0 <- xst0[k + 2]
