/* 
 * This file is part of biomcmc-lib, , a low-level library for phylogenomic analysis. 
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file distance_generator.h 
 *  \brief distance calculation between generic objects,without generating full matrix beforehand 
 */

#ifndef _biomcmc_distance_generator_h_
#define _biomcmc_distance_generator_h_

#include "distance_matrix.h" // idea is to replace/extend distance_matrix 
// TODO: how to extend to subsamples s.t. calling function doenst know but whole matrix is updated
typedef struct distance_generator_struct* distance_generator;

struct distance_generator_struct
{
  int n_samples, n_distances; // how many elements (samples) in matrix, and how many distances the function calculates at once
  int which_distance;  // which of the n_distances is being currently used
  double **dist; // distance, allowing for negative values
  bool *cached;  // if pair has been calculated or not (we assume all distances for this pair are calculated together)
  void *data;    // extra data (original features, sequences, etc. used by the distance_function() )
  void (*distance_function) (void*, int, int, double*); // defined elsewhere, receives data, i, and j, returns double[]
  int ref_counter;
};

distance_generator new_distance_generator (int n_samples, int n_distances);
void del_distance_generator (distance_generator d);
double distance_generator_get_at_distance (distance_generator d, int i, int j, int which_distance);
double distance_generator_get (distance_generator d, int i, int j);
/*! \brief defines distance calculation function wrapper, and all extra data needed by wrapper; no check is done here, but
 * wrapper should return at least as many distances sd n_distances (wrapper functions can check) */
void distance_generator_set_function_data (distance_generator d, void (*lowlevel_dist_funct)(void*, int, int, double*), void *extra_data);
/*! \brief distance wrapper may return several distances, but only one is returned by get(); this sets which
 * one (should be called before e.g. clustering) */
void distance_generator_set_which_distance (distance_generator d, int which_distance);
void distance_generator_reset (distance_generator d);

#endif
