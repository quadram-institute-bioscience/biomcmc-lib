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

#include "distance_matrix.h"

distance_matrix
new_distance_matrix (int nseqs)
{
  distance_matrix dist;
  int i, j;

  dist = (distance_matrix) biomcmc_malloc (sizeof (struct distance_matrix_struct));
  dist->ref_counter = 1;
  dist->size = nseqs;
  dist->d = (double**) biomcmc_malloc (nseqs * sizeof (double*));
  for (i=0; i < nseqs; i++) {
    dist->d[i] = (double*) biomcmc_malloc (nseqs * sizeof (double));
    /* value initialization (innocuous in sequence pairwise distances, helpful in species reconstruction) */
    dist->d[i][i] = 0.; /* diagonal elements (unused, usually) */
    for (j = 0;    j < i;     j++) dist->d[i][j] = -1.e35; /* lower triangular */
    for (j =i + 1; j < nseqs; j++) dist->d[i][j] = 1.e35;  /* upper triangular */
  }
  for (i=0; i < 20; i++) dist->freq[i] = 0.;
  dist->mean_JC_dist = dist->mean_K2P_dist = dist->mean_R = dist->var_K2P_dist = dist->var_R = 0.;

  dist->fromroot = NULL; /* allocated only if patristic distances for a topology are wanted) */
  dist->idx = dist->i_l = dist->i_r = NULL; /* idx are leaves as they appear in postorder */

  return dist;
}

void
zero_lower_distance_matrix (distance_matrix dist)
{
  int i, j;
  for (i = 1; i < dist->size; i++) for (j = 0; j < i; j++) {
    dist->d[i][j] = 0.; /* lower triangular can be used for mean values (in gene/sptree distances) */
    dist->d[j][i] = 1.e35; /* upper triangular is usually for minimum values */
  }
}

void
transpose_distance_matrix (distance_matrix dist)
{
  int i, j;
  double tmpdist;
  for (i = 1; i < dist->size; i++) for (j = 0; j < i; j++) { 
    tmpdist = dist->d[i][j];
    dist->d[i][j] = dist->d[j][i];
    dist->d[j][i] = tmpdist;
  }
}

void
del_distance_matrix (distance_matrix dist)
{
  int i;
  if (!dist) return;
  if (--dist->ref_counter) return;
  if (dist->d) {
    for (i = dist->size-1; i >= 0; i--) if (dist->d[i]) free (dist->d[i]);
    free (dist->d);
  }
  if (dist->fromroot) free (dist->fromroot);
  if (dist->idx)      free (dist->idx); /* the others (i_l and i_r) are pointers to idx elements */
  free (dist);
}

/* species-based distance matrix functions */

spdist_matrix
new_spdist_matrix (int n_species)
{
  spdist_matrix dist;
  int i, n_pairs = (n_species - 1)*(n_species)/2;

  dist = (spdist_matrix) biomcmc_malloc (sizeof (struct spdist_matrix_struct));
  dist->ref_counter = 1;
  dist->size = n_species;
  dist->n_missing = n_pairs; // number of missing comparisons
  dist->mean  = (double*) biomcmc_malloc (n_pairs * sizeof (double));
  dist->min   = (double*) biomcmc_malloc (n_pairs * sizeof (double));
  dist->count = (int*)    biomcmc_malloc (n_pairs * sizeof (int));
  dist->species_present = (bool*) biomcmc_malloc (n_species * sizeof (bool));
  for (i=0; i < n_pairs; i++) {
    dist->mean[i] = 0.; dist->min[i] = 1.e35; dist->count[i] = 0;
  }
  for (i=0; i < n_species; i++) dist->species_present[i] = false;
  return dist;
}

void
zero_all_spdist_matrix (spdist_matrix dist)
{ /* zero both mean and min, since assumes this is across loci (and only means are accounted) */
  int i, n_pairs = (dist->size - 1)*(dist->size)/2;
  dist->n_missing = n_pairs; // number of missing comparisons
  for (i=0; i < n_pairs; i++) {
    dist->mean[i] = 0.; dist->min[i] = 0.; dist->count[i] = 0;
  }
  for (i=0; i < dist->size; i++) dist->species_present[i] = false;
  return;
}

void
finalise_spdist_matrix (spdist_matrix dist)
{ /* calculate averages across loci over within-locus means and within-locus mins */
  int i, n_pairs = (dist->size - 1)*(dist->size)/2;
  double max_mean = -1.e35, max_min = -1.e35;
  for (i=0; i < n_pairs; i++) if (dist->count[i]) {
    dist->n_missing--; // one less missing pairwise comparison
    dist->mean[i] /= (double)(dist->count[i]);
    dist->min[i]  /= (double)(dist->count[i]); /* reminder: min is within locus, b/c across loci is always average */
    if (max_mean < dist->mean[i]) max_mean = dist->mean[i];
    if (max_min < dist->min[i]) max_min = dist->min[i];
    dist->count[i] = 1; /* we don't need to know it anymore, but may use when resampling/averaring several matrices */
  }
  for (i=0; i < n_pairs; i++) if (dist->count[i]) { // rescale all values to one (except missing values, which will be 1.00001)
    dist->mean[i] /= max_mean;
    dist->min[i] /= max_min;
  }
  if (dist->n_missing) for (i=0; i < n_pairs; i++) if (! dist->count[i]) {
    dist->mean[i] = 1.00001;
    dist->min[i]  = 1.00001;
  }
}

void
complete_missing_spdist_from_global_spdist (spdist_matrix local, spdist_matrix global)
{
  int i, n_pairs = (local->size - 1)*(local->size)/2;
  for (i = 0; i < n_pairs; i++) if (! local->count[i]) {
    local->mean[i] = global->mean[i];
    local->min[i] = global->min[i];
    local->count[i] = global->count[i];  // could also be zero BTW
    if (local->count[i]) local->n_missing--;
  }
  for (i=0; i < local->size; i++) if (!local->species_present[i]) local->species_present[i] = global->species_present[i];
}

void
copy_spdist_matrix_to_distance_matrix_upper (spdist_matrix spd, distance_matrix dist, bool use_means)
{ /* UPGMA< bioNJ work with upper diagonal of a square distance matrix */
// upper diagonal --> i < j in d[i][j]; index in 1d vector--> j(j-1)/2 + i 
  int i,j;
  double *sp_dist = spd->min; // default is to use min dists

  if (spd->size != dist->size) biomcmc_error ("distance matrix for NJ and species-based spdist_matrix have different sizes\n");
  if (use_means) sp_dist = spd->mean; /* alternative to use one or another would be to fill both lower and upper of dist (but a biy more expensive later to transpose) */
  for (j = 1; j < spd->size; j++) for (i = 0; i < j; i++) dist->d[i][j] = sp_dist[ (j * (j-1) / 2 + i) ];
  return;
}

// TODO: create update_spdist() to get missing vals from outside matrix; also create average_spdist(dist1, dist2)

void
del_spdist_matrix (spdist_matrix dist)
{
  if (!dist) return;
  if (--dist->ref_counter) return;
  if (dist->mean) free (dist->mean); 
  if (dist->min) free (dist->min); 
  if (dist->count) free (dist->count); 
  if (dist->species_present) free (dist->species_present);
  free (dist);
}

