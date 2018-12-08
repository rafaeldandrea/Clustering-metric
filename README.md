# KmeansGap
Determine whether a species assemblage is clustered, plot gap curves

# Usage
KmeansGap=function(
					dat,
					nozeros=FALSE,
					multiD=FALSE,
					nullmodel='shuffle',
					numnulls=100,
					mink=1,
					maxk=NULL,
					nstartingpoints=100,
					weighting='tibshirani',
					shortcut=FALSE,
					internal.peaks=FALSE,
					plot=FALSE,
					plotquant90=TRUE,
					bands=TRUE,
					verbose=TRUE
				)

# Arguments
dat 				: data.frame	Should have an abundance column named N and either a trait column named trait 
									or multiple trait columns for higher-dimensional test.
									
									
nozeros 			: log 			If TRUE, discards species with zero abundance from the null communities.

multiD 			: log 			If FALSE, runs analysis on a single trait axis, uses the column from dat named "trait" (returns error if no such column exists).
										If TRUE, runs analysis using all columns of dat other than N (multitrait analysis).

nullmodel 		: chr 			If 'shuffle', observed abundances are permutated across observed traits. 
                   					If 'draw', null traits are drawn de novo from U(0,1), and observed abundances permutated across them. 

numnulls 		: int 			Number of null communities to test against.

mink 			: int 			Minimum number of clusters to search for.

maxk 			: int 			Maximum number of clusters to search for.

nstartingpoints 	: int			Fed to argument "centers" of function kmeans(). Number of different random starting points for the clusters.

weighting 		: chr 			Type of weight applied to within-cluster dispersion. 
										If weighting = 'tibshirani', uses dispersion as defined in Tibshirani et al 2001, J R Stat Soc B 63, part 2, pp 411-423.
										If weighting = 'yan', uses correction proposed in Yan & Ye 2007, Biometrics 63, 1031-1037.

shortcut 		: log 			If TRUE, initial centroids are evenly distributed. If FALSE, initial centroids are randomly distributed, 
										and metric repeats analysis for a total of nstartingpoints times. 
										TRUE runs faster, but increases the risk of missing the global optimum.

internal.peaks 	: log			If TRUE, the metric ignores the option of a single cluster, looking instead for internal peaks in the gap curve.
									Use when looking specifically for substructure inside a single community-wide cluster, 
										as e.g. caused by environmental filters for intermediate traits. 

plot 			: log 			If TRUE, plots results as a a gap curve.

bands			: log			If TRUE, adds 95th quantile of the distribution of gap index across null communities for each number of clusters.

plotquant90 		: log 			If TRUE, adds 90th quantile of the null communities.

verbose 			: log 			If TRUE, prints dots on console indicating which number of clusters is currently being tested.

# Value
data				: data.frame	Rows show values corresponding to each number of clusters tested. Columns are as follows
			k			: num number of clusters tested (all integers between mink and maxk)
										
                        gap			: num gap index = difference in log dispersal between observed community and mean of null communities
										
                        Egap		: num mean gap index across null communities 
				 						
                        sdgap		: num standard deviation of the gap index across null communities 
										
                        nullquant95	: num 95th quantile of the gap index across null communities
										
                        nullquant90	: num 90th quantile of the gap index across null communities
										
                        logWk		: num natural logarithm of the within-cluster dispersal returned from kmeans()
										
                        ElogWk		: num mean of the values of logWk across the null communities

nullmodel		: chr 			Either 'draw' or 'shuffle'

khat				: int 			Number of clusters estimated for the observed community = number of clusters at which the gap index was maximal

maxgap			: num 			Gap statistic = maximum value of the gap index across all number of clusters tested in the observed community

maxnullquant95	: num 			95th quantile of the gap statistics of the null communities

maxnullquant90	: num			90th quantile of the gap statistics of the null communities

mink				: num 			Minimum number of clusters tested

maxk				: num 			Maximum number of clusters tested

z.score			: num 			Defined as (maxgap-mean(maxnullgap))/sd(maxnullgap)

p.value			: num 			Fraction of null communities with a higher gap statistic than the observed community.

weighting		: chr 			Prints the input value of weighting. See ARGUMENTS.

dispersal		: chr 			Prints the corresponding dispersion weighting. Let D be the within-cluster sum of squared distances:
										D/2n is used in Tibshirani et al. 2001. 
									D/(2n(n-1)) is used in Yan & Ye 2007. 

# Example
set.seed(0)

com = sample(seq(0, 1, l = 100), size = 2000, replace = TRUE, prob = 1 + sin(4 * 2 * pi * trait))

dat = plyr::count(com); names(dat) = c('trait', 'N')

plot(dat, t = 'h', las = 1)

gap = KmeansGap(dat)


# References 
Gap statistic method: R. Tibshirani, G. Walther, and T. Hastie. (2001) Estimating the number of clusters in a data set via the gap statistic. J. R. Statist. Soc. B 63, Part 2, pp. 411-423

R. D'Andrea, M. Riolo, and A. Ostling. (2018) Generalizing clusters of similar species as a signature of coexistence under competition. biorXiv, https://www.biorxiv.org/content/early/2018/02/13/264606 
