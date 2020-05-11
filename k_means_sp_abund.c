/*****
 ** MADE MANY MORE MODIFICATIONS THAN LISTED HERE, BUT KEEPING COMMENTS SO KNOW WHAT STARTED WITH
 ** kmeans.c
 ** - a simple k-means clustering routine
 ** - returns the cluster labels of the data points in an array
 ** - here's an example
 **   extern int *k_means(double**, int, int, int, double, double**);
 **   ...
 **   int *c = k_means(data_points, num_points, dim, 20, 1e-4, 0);
 **   for (i = 0; i < num_points; i++) {
 **      printf("data point %d is in cluster %d\n", i, c[i]);
 **   }
 **   ...
 **   free(c);
 ** Parameters
 ** - array of data points (double **data)
 ** - number of data points (int n)
 ** - dimension (int m)
 **
 ** A Ostling added data_abund array holding species abundances 
 ** One sp abund for each of  the n data points (each will be species in our application)
 ** - array of abundances (int *data_abund) of length n 
 ** In principle we could make those non-integer 
 ** but then have to change the counts array below to double 
 ** Not sure if code line 
 ** c[i][j] = counts[i] ? c1[i][j] / counts[i] : c1[i][j];
 ** relies on counts being integer as I am not used to this notation
 **
 ** - desired number of clusters (int k)
 ** - error tolerance (double t)
 **   - used as the stopping criterion, i.e. when the sum of
 **     squared euclidean distance (standard error for k-means)
 **     of an iteration is within the tolerable range from that
 **     of the previous iteration, the clusters are considered
 **     "stable", and the function returns
 **   - a suggested value would be 0.0001
 ** - output address for the final centroids (double **centroids)
 **   - user must make sure the memory is properly allocated, or
 **     pass the null pointer if not interested in the centroids
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
 ** Created on 2005-04-12 by
 ** - Roger Zhang (rogerz@cs.dal.ca)
 ** Modifications
 ** -
 ** Last compiled under Linux with gcc-3
 */
 /*
 ** src: http://cs.smu.ca/~r_zhang/code/kmeans.c
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
	  printf("Hit maximum iterations and still did not converge: dispersion = %f, old_dispersion=%f\n", dispersion, old_dispersion);
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
