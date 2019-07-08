/* 
 * This file is part of biomcmc-lib, a low-level library for phylogenomic analysis.
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * biomcmc is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */
/* code inspired by https://github.com/Michael-Gkotsis/Optics */

#include "optics_clustering.h"

optics_cluster
new_optics_cluster (int n_samples)
{
  optics_cluster oc = (optics_cluster) biomcmc_malloc (sizeof (struct optics_cluster_struct));

  oc->n_samples = n_samples;
  oc->cluster = (int *) biomcmc_malloc(n_samples * sizeof (int));
  oc->order   = (int *) biomcmc_malloc (n_samples * sizeof (int));
  oc->core   = (bool *) biomcmc_malloc (n_samples * sizeof (bool)); // 0 for order, 1 for core
  oc->core_distance  = (double *) biomcmc_malloc(n_samples * sizeof (double)); // distance for each core sample
  oc->reach_distance = (double *) biomcmc_malloc(n_samples * sizeof (double)); // reacheability distance for each sample
  // optics_cluster_reset (oc); // not needed since called by run()
  return oc;
}

void
del_optics_cluster (optics_cluster oc)
{ // TODO
}

void
optics_cluster_reset (optics_cluster oc)
{
  int i;
  for(i = 0; i < oc->n_samples; i++) {
    oc->cluster[i] = -1;
    oc->order[i] = 0;
    oc->core[i] = 0;
    oc->core_distance[i] = DBL_MAX;
    oc->reach_distance[i] = DBL_MAX;
  }
}

void
optics_cluster_run (optics_cluster oc, distance_generator dg, int minPoints, double minDist, double clustDist)
{
/*
  int dim = dimensions of each sample, int n = samples 
  int minPoints;    // amount of Minimum points for a cluster to be created
  double clustDist; // clustering distance: Minimum Distance for a point to be assigned to clusters
  double minDist;   // generating distance: Minimum Distance for a point to be assigned to clusters
  float *X =(float *)calloc(n*dim, sizeof(float *)) <- replaced by distance_generator()
*/
  int i,j,d;           // i counter for n, j counter for Neighborhood Elements, d counter for dimensions
  int cluster = 0;     // Cluster index
  int h = 0, k;             // h counter for Neighborhood Distance Calculation, k counter for files
  char c;
  int choise = 0;
  double min;
  int location, size, e  = 0;

  if (oc->n_samples != dg->n_samples) biomcmc_error ("sample sizes differ between OPTICS structure and distance_generator()");
  if (minPoints < 2) minPoints = 2;

  visited  = (bool*) biomcmc_malloc(oc->n_samples * sizeof (bool));
  seed     = (bool*) biomcmc_malloc(oc->n_samples * sizeof (bool));
  n_belong = (bool*) biomcmc_malloc(oc->n_samples * sizeof (bool));
  belong   = (bool*) biomcmc_malloc(oc->n_samples * sizeof (bool));
  num_pts  = (int *) biomcmc_malloc(oc->n_samples * sizeof (int)); // num_pts
  distance = (double *) biomcmc_malloc(oc->n_samples * sizeof (double));
  distmat  = (double *) biomcmc_malloc(oc->n_samples * oc->n_samples * sizeof (double));
  tmp_reach_d = (double *) biomcmc_malloc(oc->n_samples * sizeof (double)); // tmp_reach_d
  ord_reach_d = (double *) biomcmc_malloc(oc->n_samples * sizeof (double)); // ord_reach_d

  for(i = 0; i < oc->n_samples; i++) {
    num_pts[i] = 0; seed[i] = visited[i] = false; 
  }

  for(i = oc->n_samples - 1; i >= 0; i--) {
    if(!visited[i]) {
      num_pts[i] = 0;
      visited[i] = true;
      oc->order[h] = i;
      ord_reach_d[h] = DBL_MAX;
      for(j = oc->n_samples - 1; j >= 0; j--) {
        seed[j] =  belong[j] = false;
        distance[j] = 0.; 
        if(j != i) {
          distance[j] = distance_generator_get (dg, i, j); 
          if(distance[j] <= minDist) { belong[j] = true; num_pts[i]++;  }
        }
      }
      num_pts[i]++;
      if(num_pts[i] >= minPoints) {
        e = 0;
        for(j = oc->n_samples - 1; j >= 0; j--) if (belong[j]) {
          oc->reach_distance[j] = distance_generator_get (dg, i, j); 
          tmp_reach_d[e] = oc->reach_distance[j];
          e++;
          if(!visited[j]) seed[j] = true;
        }
        qSort(tmp_reach_d,e); // TODO 
        oc->core_distance[h] = tmp_reach_d[minPoints - 1];
        h++;
        do {
          size = 0;
          min = DBL_MAX;
          for(j = oc->n_samples - 1; j >=0; j--)  if(seed[j]) if (oc->reach_distance[j] < min) {
            min = oc->reach_distance[j];
            location = j;
          }  // STOPHERE
          Seed[location] = 0;
          visited[location] = 1;
          oc->order[h] = location;
          ord_reach_d[h] = oc->reach_distance[location];
          num_pts[location] = 0;
          h++;
          for(k = n; k--;) {
            nBelong[k] = 0; // Mark temp 0 for current Point, used for reset
            distance[k] = 0; //Reset the distance before calculating for current point
            if(k != location) {
              //Calculating and storing distances of Xi element with the rest of the points
              for(d = dim; d--;)     distance[k] += (X[k*dim + d] - X[location*dim + d])*(X[k*dim + d] - X[location*dim + d]);
              distance[k] = sqrt(distance[k]);
              if(distance[k] <= minDist) {   //Calculating the size of the Neighborhood
                nBelong[k] = 1;     //If Current elements belongs to the Neighborhood set temp 1 and Increase size
                num_pts[location]++;
              }
            }
          }
          num_pts[location]++;

          if(num_pts[location] >= minPoints) {
            e = 0;
            for(k = n; k--;) if(nBelong[k] != 0) {
              if(Seed[k] != 0) {
                double temp = oc->reach_distance[k];
                oc->reach_distance[k] = 0;
                //Calculating and storing distances of Xi element with the rest of the points
                for(d = dim; d--;) oc->reach_distance[k] += (X[k*dim + d] - X[location*dim + d])*(X[k*dim + d] - X[location*dim + d]);
                oc->reach_distance[k] = sqrt(oc->reach_distance[k]);
                if(temp < oc->reach_distance[k]) oc->reach_distance[k] = temp;
                tmp_reach_d[e] = oc->reach_distance[k];
                e++;
              } 
              else {
                oc->reach_distance[k] = 0;
                //Calculating and storing distances of Xi element with the rest of the points
                for(d = dim; d--;)  oc->reach_distance[k] += (X[k*dim + d] - X[location*dim + d])*(X[k*dim + d] - X[location*dim + d]);
                oc->reach_distance[k] = sqrt(oc->reach_distance[k]);
                tmp_reach_d[e] = oc->reach_distance[k];
                e++;
                if(visited[k] != 1) Seed[k] = 1;
              }
            }
            qSort(tmp_reach_d,e);
            coreDistance[h] = tmp_reach_d[minPoints - 1];
          }
          for(j = n; j--;) if(Seed[j] == 1) size++;
        } while (size != 0);
      }
      else {
        h++;
        coreDistance[h] = DBL_MAX;
      }
    }
  }
  cluster = 0;
  for(j = h; j--; ) {
    if (ord_reach_d[j] > clustDist) {
      if (coreDistance[j] <= clustDist) {
        cluster++;
        Cluster[j] = cluster;
        Noise[j] = 0;
      } 
      else  Noise[j] = 1;
    } 
    else {
      Cluster[j] = cluster;
      Noise[j] = 0;
    }
  }
  /* ---------------------------------END-------------------------------------- */
  // for(i = h; i--;) if(Noise[i] == 1) CounterNoise++;
  // printf("Clusters : %d \n",cluster );
  // printf("CounterNoise: %d\n",CounterNoise );
  // for(j = h; j--;) if(Noise[j] == 1) for(d = dim; d--;) fprintf(NoiseFile, "%lf ",OrderList[j*dim + d] );
  // for(j = h; j--;) fprintf(ReachabilityDistFile, "%lf \n",ord_reach_d[j]);
  // for(j = h; j--;) fprintf(CoreFile, "%lf \n",coreDistance[j]);
  // for(j = h; j--;) for(d = dim; d--;) fprintf(OrderFile, "%lf ",OrderList[j*dim + d]);

  free(X);
  free(Cluster);
  free(visited);
  free(Noise);
  free(nBelong);
  free(belong);
  free(num_pts);
  free(distance);
  free(Core);
  free(distance2);
  free(ord_reach_d);
  free(coreDistance);
  free(tmp_reach_d);
  free(oc->reach_distance);
  free(Seed);
  return 0;
}

double **transformPositive(int n,int dim, double**X)
{
  int i;
  int d;
  for(i = 0; i < n; i++) for(d = 0; d < dim; d++) X[i][d] = fabs(X[i][d]);
return X;
}

void qS(double *oc->reach_distance,int low, int high)
{
  int pi;
  if(low < high) {
    pi = partition(oc->reach_distance,low,high);
    qS(oc->reach_distance,low, pi - 1);
    qS(oc->reach_distance,pi + 1, high);
  }
}

int partition (double *oc->reach_distance,int low,int high)
{
  double pivot;
 int i;
 int j;
  pivot = oc->reach_distance[high];

  i = (low - 1);
  for(j = low; j <= high - 1; j++) if (oc->reach_distance[j] <= pivot) {
    i++;
    double temp = oc->reach_distance[i];
    oc->reach_distance[i] = oc->reach_distance[j];
    oc->reach_distance[j] = temp;
  }
  double temp = oc->reach_distance[i + 1];
  oc->reach_distance[i + 1] = oc->reach_distance[high];
  oc->reach_distance[high] = temp;
  return(i + 1);
}


#include <stdio.h>
#include <stdlib.h>
#include "FileHandling.h"
#include <math.h>
#include <string.h>
#include <time.h>

// https://github.com/Michael-Gkotsis/Local_Outlier_Factor
int main(int argc, char *argv[])
{

  int dim = 0;          // Dimensions of Elements
  int n = 0;            // n Elements
  int i,j,d;            // i counter for n, j counter for Neighborhood Elements, d counter for dimensions
  int kPoints;          // Minimum points for a cluster to be created
  int h,k;              // h counter for every other Element, k counter for files
  char filename[40];    // The given file
  clock_t start, end;   // Variables for counting execution time

  double reachDist = 0;
  double max;
  int choise;// Choise holder for switch
  char c;

  // Reading Dataset
  FILE* Dataset;

  printf("\n Give the DataSet file:");
  scanf("%s", filename);

  Dataset = fopen(filename, "r");
  if (!Dataset)
   {
    printf("\n There is something wrong with the Dataset file! \n\n");
    return -1;
   }

  dim = getColumns(Dataset); //Getting the dimensions of each Element of the Dataset
  rewind(Dataset);


  n = getRows(Dataset);       // Getting the Elements of the Dataset
  rewind(Dataset);

  printf("\n Elements:%d \n", n-1);
  printf("\n Dimensions:%d \n", dim);
  n--;
  printf("Give the amount of K Points: ");   // kPoints
  scanf("%d",&kPoints );

  // All the necessary memory allocation
  float *X;  // Array of Elements
  X = (float*)calloc(n*dim, sizeof(float));
  float *kDistance;   // Array for holding k-Distances for the First Core point of a cluster
  kDistance = (float*)calloc(n, sizeof(float));
  float *distance;   // Array for holding Distances for each ELement with each other Element
  distance = (float*)calloc(n*n, sizeof(float));
  unsigned int *Neighborhood; //Array for holding which element belong to which Neighborhood
  Neighborhood =(int*) calloc(n*n, sizeof(int));
  float *lrd; //Array for holding  Local Reachability Density for each element
  lrd =(float *)calloc(n, sizeof(float));
  unsigned int *NeighborhoodSize; //Array for holding the size of each Neighborhood
  NeighborhoodSize = (int*)calloc(n,sizeof(int));
  float *reachDistSum; //Array for holding the sum of reachability Distances for each element
  reachDistSum = (float *)calloc(n,sizeof(float));
  float *LOF;//Array for holding LOF value for each Element
  LOF = (float *)calloc(n,sizeof(float));
  float *NeighborhoodLrdSum; //Array for holding the sum of each Neighborhood lrd
  NeighborhoodLrdSum =(float *) calloc(n,sizeof(float));
  float *OrderedList;
  OrderedList =(float*) calloc(n*dim,sizeof(float));
  float *tmp2;
  tmp2 =(float*) calloc(n*dim,sizeof(float*));
  float *tempDistance;   // Array for holding Distances for each ELement with each other Element
  tempDistance =(float*) calloc(n, sizeof(float));

  for(i = n; i--;)
   {
    kDistance[i] = 0;
    NeighborhoodLrdSum[i] = 0;
    reachDistSum[i] = 0;
    NeighborhoodSize[i] = 0;
    lrd[i] = 0;
   }
  // Passing elements to Array X[n][dim]
  X = getData(Dataset,n,dim,X);

  for(i = n; i--;) for(d = dim; d--;) OrderedList[i*dim + d] = X[i*dim + d];
  fclose(Dataset);
  start = clock();
  /* ---------------------------------LOF-------------------------------------- */
  //STEP 1  Finding the k-Distance of each element in the dataset
  for(i = n; i--;) {
    distance[i*n + i] =  9999;
    tempDistance[i] = 99999;
    for(h = n; h--;) if(h != i) {
      distance[i*n + h] =  0;
      /* Calculating the distance of each element with i element
         and then we store it to an array so that we can order it later,
         we set the distance of i to max value so that we can avoid it at
         the ordering */
      for(d = dim; d--;) distance[i*n + h] += (X[h*dim + d] - X[i*dim + d])*(X[h*dim + d] - X[i*dim + d]);
      distance[i*n + h] = sqrt(distance[i*n +h]);
      tempDistance[h] = distance[i*n + h];
    }
    /* Ordering the array holding the distances from min to max,
       (distance[i] will be located at the distance[max]) */
    qSort(tempDistance,n);
    /* Using the ordered array to extract the k-distance of i element,
       DIST_k(i) */
    kDistance[i] = tempDistance[kPoints - 1];
  }
  //STEP 2 Finding the Neighborhood of each element
  for(i = n; i--;) {
    Neighborhood[i*n + i] = 0;
    for(h = n; h--;) if(h != i) {
      /* Calculate distance of each element with i element then store it to an array so we can determine if that element belongs to the Neighborhood of i element */
      Neighborhood[i*n + h] = 0;
      /* Comparing distance of each potential Neighbor with the k-distance of i element; if less, set Neighborhood[master][currentElement] to 1 (it belongs to i element's Neighborhood) */
      if(distance[i*n + h] <= kDistance[i])  Neighborhood[i*n + h] = 1;
    }
  }
  //STEP 3 Getting the Neighborhood size of each element
  for(i = n; i--;) {
    unsigned int counter = 0;
    for(h = n; h--;)
      if(Neighborhood[i*n + h] != 0)  counter++;
    NeighborhoodSize[i] = counter;
  }
  // step 4: Finding the Sum of reachDist_k(h <- i)
  for(i = n; i--;) for(h = n; h--;) if(Neighborhood[i*n + h] != 0) { /*Calculating the reach dist of i in respect to h so we can compare it with the k-distance of h */
    reachDist = 0;
    for(d = dim; d--;) reachDist += (X[i*dim + d] - X[h*dim + d])*(X[i*dim + d] - X[h*dim + d]);
    reachDist = sqrt(reachDist);
    /*Comparing kDistance of h w/ reach dist of i in respect to h, we choose the Max between them, max(kDistance(h),reachDist_k(h <- i))*/
    if(kDistance[h] > reachDist) max = kDistance[h];
    else max = reachDist;
    reachDistSum[i] += max;  //  Adding the max to the total Sum
  }
  /* step 4: Calculating the Local Reachability Density(lrd) of each element lrd_k(i) = NeighborhoodSize(i)/reachDistSum(i) */
  for(i = n; i--;)  lrd[i] = NeighborhoodSize[i]/(reachDistSum[i]);
  /* step 5: Calculating the sum of a Neighborhood element lrd divided by the master's element lrd for each i. Sum(lrd[h]/lrd[i]) for each Neighborhood element */
  for(i = n; i--;) for(h = n; h--;) if(Neighborhood[i*n + h] != 0) NeighborhoodLrdSum[i] += lrd[h]/lrd[i];
  //STEP 6  //Calculating the LOF for each element
  for(i = n; i--;) LOF[i] = (NeighborhoodLrdSum[i]/NeighborhoodSize[i]);
  //STEP 7  Ordering the data according to LOF from max to min
  for(i = n; i--; ) for(h = n; h--;) if(LOF[i] >= LOF[h]) {
    double tmp = LOF[i];
    LOF[i] = LOF[h];
    LOF[h] = tmp;
    for(d = dim; d--;) {
      tmp2[i*dim + d] = OrderedList[i*dim + d];
      OrderedList[i*dim + d] = OrderedList[h*dim + d];
      OrderedList[h*dim + d] = tmp2[i*dim + d];
    }
  }
  /* ---------------------------------END-------------------------------------- */
  //for(i = n; i--;) {
  //  fprintf(OrderFile, "LOF : %lf --> :",LOF[i]);
  //  for(d = dim; d--;) fprintf(OrderFile,"%lf ",OrderedList[i*dim + d] );
  //  fprintf(OrderFile, "\n");
  // }

  fclose(OrderFile);

  free(X);
  free(kDistance);
  free(distance);
  free(Neighborhood);
  free(lrd);
  free(NeighborhoodSize);
  free(reachDistSum);
  free(NeighborhoodLrdSum);
  free(LOF);
  free(tmp2);
  free(OrderedList);
  free(tempDistance);
  return 0;

}
