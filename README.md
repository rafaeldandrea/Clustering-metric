# KmeansGap
Determine whether a species assemblage is clustered, plot gap curves

# Usage
KmeansGap(dat, nozeros = FALSE, multiD = FALSE, numnulls = 100, mink = 1, maxk = NULL, nstartingpoints = 100, weighting = 0, plot = FALSE, plotquant90 = TRUE, verbose = TRUE)

# Arguments
dat: data frame with an abundance column named N and either a trait column named trait or multiple trait columns for higher-dimensional tests.

nozeros: logical. Discard species with zero abundance from the null communities?

multiD: logical. Perform analysis using all columns of dat other than N (multitrait analysis)? If FALSE, use just the one named "trait".

numnulls: integer. Number of null communities to test against.

mink: integer. Minimum number of clusters to search for.

maxk: integer. Maximum number of clusters to search for.

nstartingpoints: integer. fed to argument "centers" of function kmeans(). Number of different random starting points for the clusters.

weighting: integer. possible values are 0, 1, 2. Type of weight applied to within-cluster dispersal. 0 corresponds to no weighting. For options 1 and 2, see Tibshirani et al. 2001 and Yan & Ye 2007, respectively.

plot: logical. Plot results as a a gap curve?

plotquant90: logical. If plotting results, plot 90th quantile of the null communities?

verbose: logical. Print dots on console indicating which number of clusters is currently being tested?


# Value
KmeansGap returns a list with a data frame "data" and the summary statistics

data: data frame. Rows show values corresponding to each number of clusters tested. Columns are as follows

$k: number of clusters tested (all integers between mink and maxk)

$gap: gap index = difference in log dispersal between observed community and mean of null communities

$Egap: mean gap index across null communities 

$sdgap: standard deviation of the gap index across null communities 

$nullquant: 95th quantile of the gap index across null communities

$nullquant90: 90th quantile of the gap index across null communities

$kmaxnullquant: 95th quantile of gap index, taken only among those null communities whose max gap occurred at k clusters

$kmaxnullquant90: 90th quantile of gap index, taken only among those null communities whose max gap occurred at k clusters

$logWk: log of the within-cluster dispersal returned from kmeans()

$ElogWk: mean of the values of logWk across the null communities


khat: integer. Number of clusters estimated for the observed community = number of clusters at which the gap index was maximal

maxgap: real. Gap statistic = maximum value of the gap index across all number of clusters tested in the observed community

maxnullquant: real. 95th quantile of the gap statistics of the null communities

maxnullquant90: real. 90th quantile of the gap statistics of the null communities

mink: integer. Minimum number of clusters tested

maxk: integer. Maximum number of clusters tested

z.score: real. Defined as (maxgap-mean(maxnullgap))/sd(maxnullgap)

p.value: real. Fraction of null communities with a higher gap statistic than the observed community.

dispersal: character indicating which weighting was used. 'D' corresponds to no weighting, 'D/2n' is the weighting in Tibshirani et al. 2001, and 'D/(2n(n-1))' is the weighting in Yan & Ye 2007, where D is the within-cluster dispersal from the kmeans function.


# Example
set.seed(0)

com = sample(seq(0, 1, l = 100), size = 2000, replace = TRUE, prob = 1 + sin(4 * 2 * pi * trait))

dat = plyr::count(com); names(dat) = c('trait', 'N')

plot(dat, t = 'h', las = 1)

gap = KmeansGap(dat)


# References 
Gap statistic method: R. Tibshirani, G. Walther, and T. Hastie. (2001) Estimating the number of clusters in a data set via the gap statistic. J. R. Statist. Soc. B 63, Part 2, pp. 411-423

R. D'Andrea, M. Riolo, and A. Ostling. (2018) Competition and immigration lead to clusters of similar species, not trait separation. biorXiv, https://www.biorxiv.org/content/early/2018/02/13/264606 
