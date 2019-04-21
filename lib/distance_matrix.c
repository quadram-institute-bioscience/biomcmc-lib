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
// TODO: clann estimates missing dist pairs from dists between them and common leaf

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
    dist->mean[i] = 0.; dist->min[i] = DBL_MAX; dist->count[i] = 0;
  }
  for (i=0; i < n_species; i++) dist->species_present[i] = false;
  return dist;
}

void
zero_all_spdist_matrix (spdist_matrix dist, bool is_global)
{ /* zero both mean and min, since assumes this is across loci (and only means are accounted) */
  int i, n_pairs = (dist->size - 1)*(dist->size)/2;
  double min_value = DBL_MAX;
  if (is_global) min_value = 0.;
  dist->n_missing = n_pairs; // number of missing comparisons
  for (i=0; i < n_pairs; i++) {
    dist->mean[i] = 0.; dist->min[i] = min_value; dist->count[i] = 0;
  }
  for (i=0; i < dist->size; i++) dist->species_present[i] = false;
  return;
}

void
finalise_spdist_matrix (spdist_matrix dist)
{ /* calculate averages across loci over within-locus means and within-locus mins */
  int i, n_pairs = (dist->size - 1)*(dist->size)/2;
  double max_mean = DBL_MIN, max_min = DBL_MIN;
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
    dist->mean[i] = 1.0001;
    dist->min[i]  = 1.0001;
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
{ /* UPGMA, bioNJ work with upper diagonal of a square distance matrix */
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

void
fill_species_dists_from_gene_dists (distance_matrix spdist, distance_matrix gendist, int *sp_id, bool use_upper_gene)
{
  int i, j, k, row, col, *freq;

  freq = (int*) biomcmc_malloc (spdist->size * sizeof (int));
  for (i = 0; i < spdist->size; i++) freq[i] = 0; /* species frequency for this gene */
  for (i = 0; i < gendist->size; i++) freq[ sp_id[i] ]++; /* used to calculate mean */
  for (i = 0; i < spdist->size; i++) {
    for (j = 0; j <= i; j++)     spdist->d[i][j] = 0.; /* lower diag are mean values */
    for (;j < spdist->size; j++) spdist->d[i][j] = 1.e35; /* upper diag are minimum values */
  }
  
  for (j=1; j < gendist->size; j++) for (i=0; i < j; i++) if (sp_id[i] != sp_id[j]) {
    if (sp_id[i] < sp_id[j]) { row = sp_id[i]; col = sp_id[j]; } /* [row][col] of sptree is upper triangular for minimum */
    else                     { row = sp_id[j]; col = sp_id[i]; }
    if (!use_upper_gene) { k = i; i = j; j = k; } /* then i should be larger than j -- swap values */
    if (gendist->d[i][j] < spdist->d[row][col]) spdist->d[row][col] = gendist->d[i][j]; /* upper diag = minimum */
    spdist->d[col][row] += gendist->d[i][j]; /* lower diag = mean */
  }

  for (i = 0; i < spdist->size; i++) for (j = 0; j < i; j++) if (freq[i] && freq[j]) spdist->d[i][j] /= (double)(freq[i] * freq[j]); 

#ifdef BIOMCMC_PRINT_DEBUG
  for (i=0; i < gendist->size; i++) printf ("spdistfromgene %d\t -> %d\n", i, sp_id[i]);
  for (j=1; j < spdist->size; j++)  for (i=0; i < j; i++) 
    printf ("spdistfromgene (%d\t%d)\t%lf\n", i, j, spdist->d[i][j]); 
#endif
  free (freq);
}

void
update_species_dists_from_spdist (distance_matrix global, distance_matrix local, int *spexist)
{ 
  int i, j;
  if (global->size != local->size) biomcmc_error ("species distance matrices have different sizes within and across loci");

  for (i = 0; i < local->size; i++) for (j = 0; j < i; j++) if (spexist[i] && spexist[j]) { 
    if (global->d[j][i] > local->d[j][i]) global->d[j][i] = local->d[j][i]; /* upper triangular => minimum */
    global->d[i][j] += local->d[i][j]; /* just the sum; to have the mean we need to divide by representativity of each species across loci */
    // // guenomu receives another matrix // if (counter) { counter->d[i][j] += 1.; counter->d[j][i] += 1.; }
  }
}

int
prepare_spdistmatrix_from_gene_species_map (spdist_matrix spdist, int *sp_id, int n_sp_id)
{
  int i, number_of_species_present_in_gene = 0;
  for (i = 0; i < spdist->size; i++) spdist->species_present[i] = false; // update presence mask of species in gene 
  for (i = 0; i < n_sp_id; i++) spdist->species_present[ sp_id[i] ] = true;
  for (i = 0; i < spdist->size; i++) if (spdist->species_present[i]) number_of_species_present_in_gene++;
  return number_of_species_present_in_gene;
}

void 
fill_spdistmatrix_from_gene_dists (spdist_matrix spdist, distance_matrix gendist, int *sp_id, bool use_upper_gene)
{ // more compact than functions above (which could be eliminated in future versions)
  int i, j, i2, j2, idx, row, col, n_pairs = spdist->size*(spdist->size-1)/2;

  for (i = 0; i < n_pairs; i++) {
    spdist->mean[i] = 0;
    spdist->min[i] = 1.e35;
    spdist->count[i] = 0;
  }
  
  for (j=1; j < gendist->size; j++) for (i=0; i < j; i++) if (sp_id[i] != sp_id[j]) {
    if (sp_id[i] < sp_id[j]) { row = sp_id[i]; col = sp_id[j]; } /* make sure that row < col */
    else                     { row = sp_id[j]; col = sp_id[i]; }
    i2 = i; j2 = j;  // i2 < j2 if upper and i2 > j2 if lower diag is used
    if (!use_upper_gene) { i2 = j; j2 = i; } /* then i2 should be larger than j2 -- swap values */
    idx = col * (col-1)/2 + row; /* index in spdist */
    if (gendist->d[i2][j2] < spdist->min[idx]) spdist->min[idx] = gendist->d[i2][j2];
    spdist->mean[idx] += gendist->d[i2][j2];
    spdist->count[idx]++;
  }

  for (i = 0; i < n_pairs; i++) if (spdist->count[i]) spdist->mean[i] /= spdist->count[i];
  return;
}

void 
fill_spdistmatrix_from_gene_dist_vector (spdist_matrix spdist, double *gdist, int n_gdist, int *sp_id)
{ 
  int i, j, idx_s, idx_g, row, col, n_pairs = spdist->size*(spdist->size-1)/2;

  for (i = 0; i < n_pairs; i++) {
    spdist->mean[i] = 0;
    spdist->min[i] = DBL_MAX;
    spdist->count[i] = 0;
  }
  for (i = 0; i < spdist->size; i++) spdist->species_present[i] = false; 
  for (i = 0; i < n_gdist; i++) spdist->species_present[ sp_id[i] ] = true;

  for (j=1; j < n_gdist; j++) for (i=0; i < j; i++) if (sp_id[i] != sp_id[j]) {
    if (sp_id[i] < sp_id[j]) { row = sp_id[i]; col = sp_id[j]; } /* make sure that row < col */
    else                     { row = sp_id[j]; col = sp_id[i]; }
    idx_s = (col * (col-1))/2 + row; /* index in spdist */
    idx_g = (j * (j - 1))/2 + i; /* index in gene vector */
    if (gdist[idx_g] < spdist->min[idx_s]) spdist->min[idx_s] = gdist[idx_g];
    spdist->mean[idx_s] += gdist[idx_g];
    spdist->count[idx_s]++;
  }

  for (i = 0; i < n_pairs; i++) if (spdist->count[i]) spdist->mean[i] /= spdist->count[i];
  return;
}

void
update_spdistmatrix_from_spdistmatrix (spdist_matrix global, spdist_matrix local)
{ 
  int i, j, idx;
  if (global->size != local->size) biomcmc_error ("species spdist matrices have different sizes within and across loci");

  for (j = 1; j < local->size; j++) for (i = 0; i < j; i++) if (local->species_present[i] && local->species_present[j]) { 
    idx = (j * (j-1)) /2 + i; // index in 1D vector without diagonals (for diagonals replace -1 for +1 BTW) 
    global->mean[idx] += local->mean[idx]; // global only stores average across locals (min => within locus)
    global->min[idx] += local->min[idx];
    global->count[idx]++;
  }
  for (i = 0; i < global->size; i++) global->species_present[i] |= local->species_present[i]; // overall presence of species 
}

/* 
  How to use these functions for MRD (matrix repres with distances): < old way>
1.  distance_matrix genedist = new_distance_matrix_for_topology (gene_topol->nleaves);
2.  fill_distance_matrix_from_topology (genedist, gene_topol, NULL, true);
3.  fill_spdistmatrix_from_gene_dists (pool->this_gene_spdist, genedist, idx_gene_to_species, true);
4.  update_spdistmatrix_from_spdistmatrix (pool->d_total, pool->this_gene_spdist);
5.  finalise_spdist_matrix (pool->d_total);
6.  complete_missing_spdist_from_global_spdist (pool->d[i], pool->d_total);
7.  copy_spdist_matrix_to_distance_matrix_upper (dist, pool->square_matrix, use_within_gf_means);

  <new way> :  replaces 2D distance_matrix by 1D spdistmatrix, and uses branch lengths + internodal distances
1: not needed; 2: spdistmatrix; 3: spdist to spdist
*/

