/* This code is kept available the following github repository:
https://github.com/rafaeldandrea/Clustering-metric
*/

/* This code was used to carry out the kmeans algorithm on the output of stochastic niche models in the following manuscript:

D’Andrea, R**, M Riolo, and A Ostling (2019) Generalizing clusters of similar species as a signature of coexistence under competition. PLoS Computational Biology 15:e1006699. DOI: 10.1371/journal.pcbi.1006688
 */

/* Slightly modified versions of this code (to deal with two dimensional data, and rescaling of trait axes) were also used for the analyses in 

D’Andrea, R**, M Riolo, and A Ostling (In Press) Counting niches: Abundance‐by‐trait patterns reveal niche partitioning in a Neotropical forest. Ecology DOI: 10.1002/ecy.3019 

Those modified versions are available upon request. Email aostling@umich.edu
 */

/*****
 **
 **
 ** A Ostling created this code in 2018, using as inspiration code created on 2005-04-12 by
 ** - Roger Zhang (rogerz@cs.dal.ca)
 ** which thre seems to be a varion of acvailable at:
 ** http://labshare.cshl.edu/shares/hannonlab/www-data/DelasVives_RoeTFs_mm9/libbeato-master/beato/cluster.c
 ** 
 **The code I had stated the location http://cs.smu.ca/~r_zhang/code/kmeans.c but it no longer seems to be available there
 **
 **
 ** That code had the following References and Notes:
 ** References
 ** - J. MacQueen, "Some methods for classification and analysis
 **   of multivariate observations", Fifth Berkeley Symposium on
 **   Math Statistics and Probability, 281-297, 1967.
 ** - I.S. Dhillon and D.S. Modha, "A data-clustering algorithm
 **   on distributed memory multiprocessors",
 **   Large-Scale Parallel Data Mining, 245-260, 1999.
 ** Notes
 ** - this function is provided as is with no warranty.
 ** - the author is not responsible for any damage caused
 **   either directly or indirectly by using this function.
 ** - anybody is free to do whatever he/she wants with this
 **   function as long as this header section is preserved.
 **
 **
 */
/*A Ostling Notes*/
 /* This particular code reads in a list of species' trait values and abundances. It is set up for a single trait--the analysis is one dimensional, though it could easily be generalized as in the cluster.c code mentioned above. */

/* It is actually set up to run on a cluster and read in a series of such files and analyze them, ranging from stnum_NULL.txt to ednnum_NULL.txt. Also, if it is instructed to read in 1_NULL.txt it will also analyze 0_NULL.txt. This structure was set up so that 0_NULL.txt was the actual data file,a nd 1_NULL.txt through x_NULL.txt where the x nulls generated from the data in some way controlled in other code. In some other uses of this kmeans code, for analysis of field data, a given job would only analyze part of the nulls for the data. As you will see for the cluster scripts used for the analysis of outputs of simulations, we used 100 nulls and analyzed them all in one job.*/

/* For a given trait abundance file it carries out the kmeans algorithm to find the best arrangment of species into a specified number of clusters. It interprets each species as n points at the species' trait value for the purposes of the kmeans calculations, but manipulates those points all together for efficiency (where normally the kmeans algorithm works point by point because when it is not species' abundance data it is unlikely for two points to be at exactly the same position). In this case it looks for the best arrangement for between 1 and 20 clusters. 

   /* Note our approach takes a user specified of starting points for the cluster centers that depends on the number of clusters being searched for. These are listed in a file nstart.txt We determined these nubmers of starting points needed in advance from test runs to see when the dispersion index became saturated as the number of starting points was increased more. This did seem to depend a bit on the particular trait abundance data, but we chose values that seemed to give saturation across all of the community dynamics we ran in the simulations. We found that using a substantial number of staring points was important for cluster detection.*/

/* Note for each run of the kmeans algorithm from a given starting point it continues to look for a cluster arrangement that is better until the difference in dispersion from one iteration through the algorithm to the next is smaller than a tolerance t (which we normally set to 0.0000001), or until it has done iter_MAX iterations.*/

/* Note that because this was set up to run on a cluster as part of a large job array, it was also set to read the particular data directory to use for this particular job from a path_file.*/

/* Some notes about the dispersion metrics calculated in the code: */

/*"best_dispersion" is the Wk defined in Eq. 2 of Tibhsirani et al. 2001. Note their explanation of the formula is a bit confusing. The formula shows Wk defined as the sume over clusters of 1/2n_r Dr, where Dr is the sum of distances between all points within the rth cluster and n_r is the number of points in it. However, then it is noted that Wk is also the sum over clusters of the squared distances from the cluster center, which is the way we calculate it here. The key thing to realize to see the quivalence of these statement is that the suym over the squared distances from the cluster center is the same thing as 1/2nr Dr. Just expanding out the relevant quantites in the sum and doing some algebra will help you see it. Note that in making use of Wk and the Gap method described in Tibshirani we actually look for the maximum Gap value in the range of cluster numbers we examined, rather than the smallest number of clusters with a local peak Gap value. We found that using the peak Gap value worked better. See D'Andrea et al. 2019 and 2020 cited above for more info.*/

/* "best_dispersion_yan_ye" is the modified disperstion index suggested by Yan and Ye (2007). Specifically it is the dispersion defined in their Equation 5, but calculated by instead calcualting the sum of squared distances from all points to the center of the cluster, i.e. the equivalent of 1/2n_r Dr, and then dividing that by n_r-1 and summing over the clusters. */

/*References
Tibshirani R, Walther G, Hastie T. Estimating the number of clusters in a data set via the gap statistic.
Journal of the Royal Statistical Society: Series B (Statistical Methodology). 2001; 63:411–423. https://
doi.org/10.1111/1467-9868.00293

Yan M, Ye K. Determining the Number of Clusters Using the Weighted Gap Statistic. Biometrics 2007; 6: 1031-1037. http://www.jstor.org/stable/4541456

*/

#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

int get_rand_integ_intvl(int, int);
double get_rand_unit();
void start_time(void);
void prn_time(void);
void prn_total_time(void);

int main(int argc, char *argv[])
{
  char infile_data[2000], *infile_nstart, beg_outfile[2000], directory[2000], *path_file;
  int S,S_MAX=500;
  int max_clus_num=20;
  int nstart[max_clus_num];
  int directory_num;
  
  FILE *ifp_data, *ifp_nstart, *ifp_PathFile, *ofp_dispersion, *ofp_dispersion_yan_ye, *ofp_centroids_counts, *ofp_cluster_assignment;
  char outfile_dispersion[2000], outfile_dispersion_yan_ye[2000], outfile_centroids_counts[2000], outfile_cluster_assignment[2000];
  double *trait;
  int *abund;
  int i,j,l,m, file_num, ifile_stnum, ifile_endnum;
  char line[1024];
  const char* tok;

  double t=0.0000001;
  double *centroids;
  int k;
  int h;
  int counts[max_clus_num], **best_counts; /* size of each cluster*/
  int *labels, **best_labels; /*cluster number of each species*/
  double old_dispersion, dispersion = DBL_MAX, best_dispersion=DBL_MAX; /* sum of squared euclidean distance */
  double c[max_clus_num], c1[max_clus_num], **best_c;
  double min_distance;
  double distance;
  double old_dispersion_byclus[max_clus_num], dispersion_byclus[max_clus_num], **best_dispersion_byclus;
  double best_dispersion_yan_ye;
  int choice,already_chosen;
  int *sp_already_chosen_centroid;
  int centroid_is_diff;
  int empty_cluster;
  int num_empty_clusters;
  int iter=0, iter_MAX=100;
  
  ifile_stnum = atoi(argv[1]); /*Will start with file "stnum_NULL.txt"*/
  if (ifile_stnum==1) ifile_stnum=0; /*A kludge to get the first job in the array to also analyze the data, which is in 0_NULL.txt*/
  ifile_endnum = atoi(argv[2]); /*Will finish with file "endnum_NULL.txt"*/
  directory_num = atoi(argv[3]); /*Line in path_file to use to decide what directory to get the data file from*/
  path_file = argv[4];
  
  printf("ifile_stnum is %d\n", ifile_stnum);
  printf("ifile_endnum is %d\n", ifile_endnum);
  printf("directory_num is %d\n", directory_num);
  printf("pathfile is %s\n", path_file);
  
  printf("Will now read in the nstart values from nstart.txt.\n");
  ifp_nstart = fopen("nstart.txt", "r");
  for (k=0; k<max_clus_num; ++k) {
    fgets(line, 1024, ifp_nstart);
    tok=strtok(line,",");
    /*    printf("tok is now %s\n", tok);*/
    tok=strtok(NULL,",");
    /*    printf("tok is now %s\n", tok);*/
    sscanf(tok, "%d", &nstart[k]);
    printf("nstart for %d clusters is %d\n", k+1, nstart[k]);
  }
  fclose(ifp_nstart);

  printf("Will now read the directory path to use.\n");
  ifp_PathFile = fopen(path_file, "r");
  for (i=1; i<=directory_num; ++i) {
    fgets(line, 1024, ifp_PathFile);
    printf("current line is %s", line);
    if (i==directory_num) {
      tok=strtok(line,",");
      /*    printf("tok is now %s\n", tok);*/
      tok=strtok(NULL,",");
      /*    printf("tok is now %s\n", tok);*/
      sscanf(tok, "%s", directory);
      printf("directory to use is %s\n", directory);
    }
  }
  fclose(ifp_PathFile);

  printf("Now allocating arrays\n");
  trait = calloc(S_MAX, sizeof(double));
  abund = calloc(S_MAX, sizeof(int));
  labels = calloc(S_MAX, sizeof(int));
  best_labels = calloc(max_clus_num, sizeof(int*));
  for (k=1; k<=max_clus_num; ++k) 
	best_labels[k-1] = calloc(S_MAX, sizeof(int));
  best_counts = calloc(max_clus_num, sizeof(int*));
  best_c = calloc(max_clus_num,sizeof(double*));
  best_dispersion_byclus = calloc(max_clus_num, sizeof(double*));
  for (k=1; k<=max_clus_num; ++k) {
    best_counts[k-1] = calloc(k, sizeof(int));
    best_c[k-1] = calloc(k, sizeof(double));
    best_dispersion_byclus[k-1] = calloc(k, sizeof(double));
  }
  printf("Finished allocating \n");
  
  for (file_num=ifile_stnum; file_num<=ifile_endnum; file_num++) {
    
    printf("file_num is now %d\n", file_num);
    sprintf(infile_data, "%s/%d_NULL.txt", directory, file_num);
    sprintf(beg_outfile, "%s/output/%d_NULL", directory, file_num);
    
    printf("Will read species from %s and output to beginning outfiles: %s\n", infile_data, beg_outfile);
    

    printf("Will now read the data in file %s.\n", infile_data);
    ifp_data = fopen(infile_data, "r");
    S=0;
    while(fgets(line, 1024, ifp_data) != NULL && S<S_MAX) {
	   /*   printf("Just read line: %s\n", line);*/
	  tok=strtok(line, " ");
	  /*   printf("tok is now %s\n", tok);*/
	  sscanf(tok, "%lf", &trait[S]);
	  tok=strtok(NULL," ");
	  /*    printf("tok is now %s\n", tok);*/
	  sscanf(tok,"%d", &abund[S]);
	  /*    printf("trait=%f, abund=%d\n", trait[S],abund[S]);*/
	  if(abund[S]!=0)
		S=S+1;
    }
    fclose(ifp_data);
    printf("Number of species is %d\n", S);
    if (S==S_MAX) printf("ERROR!!! Your data may have more species than %d, fix hardcoding of array size!\n", S_MAX);

    sprintf(outfile_dispersion, "%s_dispersion.txt", beg_outfile);
    ofp_dispersion = fopen(outfile_dispersion, "w");
    fprintf(ofp_dispersion, "k,log_Wk_tibshirani,log_Wk_yanye,emptyclus,iter\n");
  
    srand(time(NULL));
    start_time();
  
    for (k=1; k<=max_clus_num; ++k) {
      best_dispersion=DBL_MAX;
      for (i=0;i<k;i++) best_dispersion_byclus[k-1][i]=0;
      num_empty_clusters=0;
      printf("about to do %d starting points for cluster number %d\n", nstart[k-1], k);

      for (l=1; l<=nstart[k-1]; ++l) {
	/*     printf("Doing starting point %d\n", l);*/
	empty_cluster=0;
      
	/*initialize cluster centers--note didn't bother to check if already chosen so may get same cluster center twice*/
	for (i=0; i<k; i++) {
	  choice=get_rand_integ_intvl(0,S-1);
	  c[i] = trait[choice];
	  /*	printf("cluster %d chose species %d at trait value %f\n", i+1, choice+1, trait[choice]);*/
	}

	iter = 0;
	dispersion=DBL_MAX;
	for (i=0;i<k;i++) dispersion_byclus[i]=0;
	do {
	  iter++;
	  centroid_is_diff=0;
	  /*	printf("Doing iteration %d\n", iter);*/
	  /* SAVE ERROR FROM LAST STEP */
	  old_dispersion = dispersion, dispersion = 0;
	  for (i=0;i<k;i++) {
	    old_dispersion_byclus[i] = dispersion_byclus[i];
	    dispersion_byclus[i] = 0;
	  }
     
	  /* CLEAR OLD COUNTS AND TEMP CENTROIDS */
	  for (i = 0; i < k; counts[i++] = 0) {
	    c1[i] =0;
	  }
	  for (h = 0; h < S; h++) {
	    /* IDENTIFY THE CLOSEST CLUSTER */
	    min_distance = DBL_MAX;
	    for (i=0; i<k; i++) {
	      distance = 0;
	      distance += pow(trait[h]-c[i],2);
	      if (distance < min_distance) {
		labels[h] = i;
		min_distance=distance;
	      }
	    }
	    /* UPDATE SIZE AND TEMP CENTROID OF THE DESTINATION CLUSTER */
	    c1[labels[h]] += trait[h]*abund[h];
	    counts[labels[h]]=counts[labels[h]] + abund[h];
	    /* UPDATE Dispersion */
	    dispersion += min_distance*abund[h];
	    dispersion_byclus[labels[h]]+= min_distance*abund[h];
	  }
	  for (i = 0; i < k; i++) { /* UPDATE ALL CENTROIDS */
	    if (counts[i] !=0 && c1[i] !=0){
	      if (c[i] != c1[i]/((double) counts[i]))centroid_is_diff=1;
	      c[i] = c1[i]/((double) counts[i]);
	    }
	    else {
	      empty_cluster = 1;
	      num_empty_clusters++;
	      /*	    printf("Cluster %d was empty\n", i);*/
	    }
	  }
	}while ((centroid_is_diff || fabs(dispersion - old_dispersion) > t) && iter < iter_MAX);
	
	/*      printf("For %d clusters, dispersion is %f, old dispersion is %f, log of dispersion is %f, iter = %d\n", k, dispersion, old_dispersion, log(dispersion), iter);*/
	/*      printf("best_dispersion was %f\n", best_dispersion);*/
	if (!(iter==iter_MAX)&&!empty_cluster && dispersion < best_dispersion) {
	  best_dispersion = dispersion;
	  for (i=0; i<k; ++i) {
	    best_c[k-1][i] = c[i];
	    best_counts[k-1][i] = counts[i];
	  }
	  for (h=0; h<S; ++h) {
	    best_labels[k-1][h] = labels[h];
	  }
	  for (i=0;i<k;++i) best_dispersion_byclus[k-1][i]=dispersion_byclus[i];
	}
	if (empty_cluster) {
	  l--;
	  /*	printf("one of the clusters was empty!\n");*/
	}
	/*      printf("it took %d iterations\n", iter);*/
	if (iter==iter_MAX) {
	  l--;
	  printf("Hit maximum iterations and still did not converge: dispersion = %f, old _dispersion=%f\n", dispersion, old_dispersion);
	}
	/*      printf("now it is %f and its log is %f\n", best_dispersion, log(best_dispersion));*/
	/*      printf("the best cluster centers are: ");
		for (i=0;i<k;++i) printf("%f,", best_c[i]);
		printf("\n");
		if (empty_cluster) printf("there was an empty cluster so setting counter back so do it again\n");
		printf("\n");*/
      }
      
      /*	      printf("---------------\n");*/
      printf("Finished with %d cluster minimization\n", k);
      prn_time();
      printf("\n\n");
      fflush(NULL);
      best_dispersion_yan_ye=0;
      for (i=0;i<k;++i) {
	best_dispersion_yan_ye+= (best_dispersion_byclus[k-1][i])/(best_counts[k-1][i]-1);
      }
      fprintf(ofp_dispersion, "%d,%f,%f,%d,%d\n",k,log(best_dispersion),log(best_dispersion_yan_ye), num_empty_clusters,iter);
    }
    sprintf(outfile_centroids_counts, "%s_centroids_counts.txt", beg_outfile);
    ofp_centroids_counts = fopen(outfile_centroids_counts, "w");
    fprintf(ofp_centroids_counts,"numclus,centroid,counts,dispersion\n");
    for (k=1;k<=max_clus_num;++k) {
      for (i=0;i<k;++i) {
		fprintf(ofp_centroids_counts, "%d,%lf,%d,%f\n",k,best_c[k-1][i], best_counts[k-1][i], best_dispersion_byclus[k-1][i]);
      }
    }
    fclose(ofp_centroids_counts); 

    sprintf(outfile_cluster_assignment, "%s_cluster_assignment.txt", beg_outfile);
    ofp_cluster_assignment = fopen(outfile_cluster_assignment, "w");
    fprintf(ofp_cluster_assignment, "trait,abund,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10,k11,k12,k13,k14,k15,k16,k17,k18,k19,k20\n");
    for (h=0;h<S;++h) {
      fprintf(ofp_cluster_assignment, "%f,%d", trait[h], abund[h]);
      for (k=1; k<=max_clus_num; ++k) 
	fprintf(ofp_cluster_assignment, ",%d", best_labels[k-1][h]+1);
      fprintf(ofp_cluster_assignment, "\n");
    }
    fclose(ofp_cluster_assignment);
    fclose(ofp_dispersion);
    prn_total_time();
  }
  
}
    
  

 /* get_rand_integ_intvl returns a random integer in the interval [x,y], that
   is, including x and y as possibilities. */
/* REMEMBER TO SEED THE RANDOM NUMBER GENERATOR when using this
   function! --srand(time(NULL)) */

int get_rand_integ_intvl(int x, int y)
{

  double toss;

  toss = get_rand_unit();
  /*  printf("toss=%f, y-x+1=%d, toss*(y-x+1)=%f, ceil(toss*(y-x+1))=%f, returning %d.\n", toss, y-x+1, toss*(y-x+1), ceil(toss*(y-x+1)), (int) ceil(toss*(y-x+1))-1 + x);*/

  /* technically should throw out zeros, because should choose a particular integer based on a region non-inclusive of the one before that, even when choosing the lowest value in the integer range */
  while (toss==0) toss=get_rand_unit();

  return (int) ceil(toss*(y-x+1))-1 + x;

  /* an incorrect modification */
  /*  return (int) ceil(toss*(y-x)) + x;*/



  /* old code that doesn't work properly: */
  /*  toss = rand();
      if (toss == 0) 
      return x;
      else
      return (ceil((toss*(y-(x-1)))/RAND_MAX) + x-1);*/

}

double get_rand_unit()
{

  double toss;

  toss = (double) rand();
  toss = toss/((double) RAND_MAX);
  
  return toss;

}

# define MAXSTRING 100

typedef struct {
  clock_t begin_clock, save_clock;
  time_t begin_time, save_time;
} time_keeper;

static time_keeper tk;

void start_time(void)
{
  tk.begin_clock = tk.save_clock = clock();
  tk.begin_time = tk.save_time = time(NULL);
}

void prn_time(void)
{
  char s1[MAXSTRING], s2[MAXSTRING];
  int field_width, n1, n2;
  double clocks_per_second = (double) CLOCKS_PER_SEC, user_time, real_time;
  
  user_time = (clock() - tk.save_clock) / clocks_per_second;
  real_time = difftime(time(NULL), tk.save_time);
  tk.save_clock = clock();
  tk.save_time = time(NULL);

  /* print the values found, and do it neatly */

  n1 = sprintf(s1, "%.1f", user_time);
  n2 = sprintf(s2, "%.1f", real_time);
  field_width = (n1 > n2) ? n1 : n2;
  printf("SINCE LAST PRN_TIME:\n");
  printf("%s%*.1f%s%f%s%f%s\n%s%*.1f%s%f%s%f%s\n", 
	 "User time: ", field_width, user_time, " seconds, ", user_time/60, " minutes, ", user_time/3600, " hours ", 
	 "Real time: ", field_width, real_time, " seconds, ", real_time/60, " minutes, ", real_time/3600, " hours ");

}

void prn_total_time(void)
{
  char s1[MAXSTRING], s2[MAXSTRING];
  int field_width, n1, n2;
  double clocks_per_second = (double) CLOCKS_PER_SEC, user_time, real_time;
  
  user_time = (clock() - tk.begin_clock) / clocks_per_second;
  real_time = difftime(time(NULL), tk.begin_time);

  /* print the values found, and do it neatly */

  n1 = sprintf(s1, "%.1f", user_time);
  n2 = sprintf(s2, "%.1f", real_time);
  field_width = (n1 > n2) ? n1 : n2;
  printf("SINCE START:\n");
  printf("%s%*.1f%s%f%s%f%s\n%s%*.1f%s%f%s%f%s\n\n", 
	 "User time: ", field_width, user_time, " seconds, ", user_time/60, " minutes, ", user_time/3600, " hours ", 
	 "Real time: ", field_width, real_time, " seconds, ", real_time/60, " minutes, ", real_time/3600, " hours ");

}
