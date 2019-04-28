/* 
 * This file is part of biomcmc-lib, a low-level library for phylogenomic analysis.
 * Copyright (C) 2019-today  Leonardo de Oliveira Martins [ leomrtns at gmail.com;  http://www.leomartins.org ]
 *
 * biomcmc is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later
 * version.

 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
 * details (file "COPYING" or http://www.gnu.org/copyleft/gpl.html).
 */

/*! \file topology_distance.h 
 *  \brief branch length operations on topologies, including patristic distances 
 */

#ifndef _biomcmc_topology_distance_h_
#define _biomcmc_topology_distance_h_

#include "topology_common.h"


/*! \brief allocate memory for a new distance_matrix that will be used on topologies */
distance_matrix new_distance_matrix_for_topology (int nleaves);

/*! \brief fill in distance_matrix with the patristic distances from topology (can be used with distinct branch length vectors to fill upper and lower diagonals */
void fill_distance_matrix_from_topology (distance_matrix dist, topology tree, double *blen, bool use_upper);

/*! \brief calculates rescaled patristic distances returning up to 6 distinct 1D vectors #dist (externally allocated) 
 * The 'tolerance' is the minimum branch length to be considered a multifurcation (length zero) */
void patristic_distances_from_topology_to_vectors (topology tree, double **dist, double *scaling, int n_dists, double tolerance);

/*! \brief similar to an Euler tour, has list of leaves below each node */
int* create_vector_with_idx_leaves_below_for_patristic (topology tree);

void estimate_topology_branch_lengths_from_distances (topology tree, double *dist);
double* new_topology_branch_lengths_from_distances (topology tree, double *dist);

#endif
