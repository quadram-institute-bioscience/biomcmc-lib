/* 
 * This file is part of biomcmc-lib, a low-level library for phylogenomic analysis.
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

#include "distance_generator.h"

distance_generator
new_distance_generator (int n_samples, int n_distances)
{
  int i, j, n_combined;
  distance_generator d = (distance_generator) biomcmc_malloc (sizeof (struct distance_generator_struct));
  if (n_distances < 1) n_distances = 1;
  d->n_distances = n_distances;
  d->n_samples = n_samples;
  n_combined = (int)(d->n_samples * (d->n_samples - 1))/2;
  d->dist = (double**) biomcmc_malloc (n_combined * sizeof (double*)); // dist[samples x samples][n_distances]
  d->cached = (bool*) biomcmc_malloc (n_combined * sizeof (bool));
  for (i=0; i < n_combined; i++) d->dist[i] = (double*)  biomcmc_malloc (d->n_distances * sizeof (double));
  for (i=0; i < n_combined; i++) {
    d->cached[i] = false; 
    for (j=0; j < d->n_distances; j++) d->dist[i][j] = 0.; 
  }
  d->data = NULL;
  d->distance_function = NULL;
  d->which_distance = 0;
  d->ref_counter = 1;
  return d;
}

void
del_distance_generator (distance_generator d)
{
  int i;
  if (!d) return;
  if (--d->ref_counter) return;
  if (d->dist) {
    for (i = (d->n_samples * (d->n_samples - 1))/2 - 1; i >=0; i--) if (d->dist[i]) free (d->dist[i]);
    free (d->dist);
  }
  if (d->cached) free (d->cached);
  free (d);
}

double
distance_generator_get (distance_generator d, int i, int j)
{
  return distance_generator_get_at_distance (d, i, j, d->which_distance);
}

double
distance_generator_get_at_distance (distance_generator d, int i, int j, int which_distance)
{
  if (i == j) return 0.;
  which_distance %= d->n_distances; // wrap around in case user gave too large which_distance
  if (j < i) { int tmp = i; i = j; j = tmp; } // upper diagonal: i<j in 2D[i][j] => 1D[j(j-1)/2 + i]
  int idx =  ((j * (j-1)) / 2 + i);
  if (! d->cached[idx]) {
    d->distance_function (d->data, i, j, d->dist[idx]); // last arg is vector where result distances will go
    d->cached[idx] = true;
  }
  return d->dist[idx][which_distance];
}

void
distance_generator_set_function_data (distance_generator d, void (*lowlevel_dist_funct)(void*, int, int, double*), void *extra_data)
{
  d->distance_function = lowlevel_dist_funct; // lowlevel_dist_funct (extra_data, i, j, *results[n_distances])
  d->data = extra_data;
}

void 
distance_generator_set_which_distance (distance_generator d, int which_distance)
{
  d->which_distance = which_distance % d->n_distances; // remainder '%' to make sure index < n_distances
}

void
distance_generator_reset (distance_generator d)
{
  int i, j, n_combined = (int)(d->n_samples * (d->n_samples - 1))/2;
  for (i=0; i < n_combined; i++) {
    d->cached[i] = false; 
    for (j=0; j < d->n_distances; j++) d->dist[i][j] = 0.; // dist can be any number actually
  }
}
