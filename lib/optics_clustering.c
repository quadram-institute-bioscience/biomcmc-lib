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
  oc->n_clusters = 0;
  oc->cluster = (int*) biomcmc_malloc (2 * n_samples * sizeof (int)); // malloc one big chunk
  oc->order   = oc->cluster + n_samples; 
  oc->core   = (bool*) biomcmc_malloc (n_samples * sizeof (bool)); // 0 for order, 1 for core
  oc->reach_distance = (double*) biomcmc_malloc (2 * n_samples * sizeof (double)); // reacheability distance for each sample
  oc->core_distance  = oc->reach_distance + n_samples;  // distance for each core sample
  return oc;
}

void
del_optics_cluster (optics_cluster oc)
{ 
  if (!oc) return;
  if (oc->reach_distance) free (oc->reach_distance);
  if (oc->core)    free (oc->core);
  if (oc->cluster) free (oc->cluster);
  free (oc);
}

void
optics_cluster_reset (optics_cluster oc)
{
  int i;
  for(i = 0; i < oc->n_samples; i++) {
    oc->cluster[i] = -1;
    oc->order[i] = 0;
    oc->core[i] = false;
    oc->core_distance[i] = DBL_MAX;
    oc->reach_distance[i] = DBL_MAX;
  }
  oc->n_clusters = 0;
}

void
optics_cluster_run (optics_cluster oc, distance_generator dg, int minPoints, double epsilon, double clustDist)
{ /* minPoints: Minimum points for a cluster to be created
     clustDist: for cluster assignment; ideal is to be decided from the reach_distance[]
     epsilon: maximum radius to consider as neighbourhood (must be large, but */
  int i, j, h = 0, k, location, size, e = 0, cluster = 0, *num_pts = NULL;
  bool *visited = NULL, *seed = NULL, *n_belong = NULL, *belong = NULL;
  double min, *distance = NULL, *tmp_reach_d = NULL, *ord_reach_d = NULL;

  if (oc->n_samples != dg->n_samples) biomcmc_error ("sample sizes differ between OPTICS structure and distance_generator()");
  if (minPoints < 2) minPoints = 2;
  optics_cluster_reset (oc);
  if (clustDist > (epsilon - 1.e-5)) clustDist = epsilon - 1.e-5; // dists higher than epsilon are not considered

  visited = (bool*) biomcmc_malloc (4 * oc->n_samples * sizeof (bool));  // alloc one big chunk
  seed     = visited +     oc->n_samples; 
  n_belong = visited + 2 * oc->n_samples;
  belong   = visited + 3 * oc->n_samples;
  num_pts  = (int*)  biomcmc_malloc (oc->n_samples * sizeof (int)); 
  distance = (double*) biomcmc_malloc (3 * oc->n_samples * sizeof (double));
  tmp_reach_d = distance +     oc->n_samples;
  ord_reach_d = distance + 2 * oc->n_samples;

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
          if(distance[j] <= epsilon) { belong[j] = true; num_pts[i]++;  }
        }
      }
      num_pts[i]++;
      if(num_pts[i] >= minPoints) {
        e = 0;
        for(j = oc->n_samples - 1; j >= 0; j--) if (belong[j]) {
          oc->reach_distance[j] = distance_generator_get (dg, i, j); 
          tmp_reach_d[e++] = oc->reach_distance[j];
          if(!visited[j]) seed[j] = true;
        }
        qsort (tmp_reach_d, e, sizeof (double), compare_double_increasing);
        oc->core_distance[h] = tmp_reach_d[minPoints - 1];
        if (oc->core_distance[h] < epsilon) oc->core[i] = true; // Leo 
        h++;
        do {
          size = 0;
          min = DBL_MAX;
          for (j = oc->n_samples - 1; j >=0; j--) if(seed[j]) if (oc->reach_distance[j] < min) {
            min = oc->reach_distance[j];
            location = j;
          }
          seed[location] = false; visited[location] = true;
          oc->order[h] = location;
          ord_reach_d[h] = oc->reach_distance[location];
          num_pts[location] = 0;
          h++;
          for(k = oc->n_samples - 1; k >= 0; k--) {
            n_belong[k] = false; 
            distance[k] = 0.;
            if(k != location) {
              distance[k] = distance_generator_get (dg, k, location); 
              if(distance[k] <= epsilon) { 
                n_belong[k] = true;
                num_pts[location]++;
              }
            }
          }
          num_pts[location]++;
          if(num_pts[location] >= minPoints) {
            e = 0;
            for(k = oc->n_samples - 1; k >=0; k--) if(n_belong[k]) {
              if(seed[k]) {
                double temp = oc->reach_distance[k];
                oc->reach_distance[k] = distance_generator_get (dg, k, location); 
                if(temp < oc->reach_distance[k]) oc->reach_distance[k] = temp;
                tmp_reach_d[e++] = oc->reach_distance[k];
              } 
              else {
                oc->reach_distance[k] = distance_generator_get (dg, k, location); 
                tmp_reach_d[e++] = oc->reach_distance[k];
                if(!visited[k]) seed[k] = true;
              }
            }
            qsort (tmp_reach_d, e, sizeof (double), compare_double_increasing);
            oc->core_distance[h] = tmp_reach_d[minPoints - 1];
            if (oc->core_distance[h] < epsilon) oc->core[location] = true; // Leo 
          }
          for(j = oc->n_samples - 1; j >= 0; j--) if (seed[j]) size++;
        } while (size != 0);
      }
      else {
        oc->core_distance[h] = DBL_MAX;
      }
    }
  }
  cluster = 0;
  for(j = h; j >= 0; j--) {
    if ((ord_reach_d[j] > clustDist) && (oc->core_distance[j] <= clustDist)) {
      cluster++;
      oc->cluster[j] = cluster;
    } 
    else oc->cluster[j] = cluster;
  }

  oc->n_clusters = cluster + 1; // cluster starts at zero, when we have _1_ cluster already

  if (distance) free (distance); // the others are just pointers to blocks of these vectors
  if (num_pts) free (num_pts);
  if (visited) free (visited);
  return;
}
