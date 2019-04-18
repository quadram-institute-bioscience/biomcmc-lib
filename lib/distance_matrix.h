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

/*! \file distance_matrix.h
 *  \brief distance matrix, that can be used in alignments and trees, and patristic-distance based species distances
 *  These functions don't know about trees/topologies: topology_common.c creates the actual patristic distances, and
 *  downstream software or genetree.h should decide how to use this information
 */

#ifndef _biomcmc_distance_matrix_h_
#define _biomcmc_distance_matrix_h_ 

#include "hashtable.h"

typedef struct distance_matrix_struct* distance_matrix;
typedef struct spdist_matrix_struct* spdist_matrix;

struct distance_matrix_struct
{
  int size;   /*! \brief number of sequences to calculate distances */
  double **d, /*! \brief pairwise distance matrix (upper) and ti/tv rate ratio (lower triangle) for K2P formula for alignments */
         mean_K2P_dist, /*! \brief average pairwise distance from K2P model */
         var_K2P_dist,  /*! \brief variance in pairwise distance from K2P model */
         mean_JC_dist,  /*! \brief average pairwise distance from JC model */
         mean_R,        /*! \brief average K2P transition/transversion ratio from pairwise distances */
         var_R,         /*! \brief variance in K2P transition/transversion ratio from pairwise distances */
         freq[20];      /*! \brief empirical equilibrium frequencies */
  double *fromroot;     /*! \brief distance from root (used to calculate distance between tree leaves) */
  int *idx, *i_l, *i_r; /*! \brief aux vectors for finding leaves spanned by subtrees on any node */
  int ref_counter;
};

struct spdist_matrix_struct
{
  int size, n_missing;
  double *mean, *min; /*! \brief mean or min distances across possibilities (within loci) */
  int *count; /*! \brief how many times this pairwise comparison appears (between or within loci) */
  bool *species_present; /*! \brief boolean marking if species is present at all in this matrix */
  int ref_counter;
};

/*! \brief creates new matrix of pairwise distances */
distance_matrix new_distance_matrix (int nseqs);
/*! \brief specially in gene/sptree distance methods (GLASS, STEAC, etc.) lower is used for means and upper for min. This function resets matrix elements */
void zero_lower_distance_matrix (distance_matrix dist);
/*! \brief invert lower and upper diagonals of matrix (since some functions like upgma expect upper, etc.) */
void transpose_distance_matrix (distance_matrix dist);
/*! \brief releases memory allocated to distance_matrix (this structure has no smart ref_counter) */
void del_distance_matrix (distance_matrix dist);

spdist_matrix new_spdist_matrix (int n_species);
void zero_all_spdist_matrix (spdist_matrix dist); /* zero both mean[] and min[] since we only look at AVERAGE across loci (never the min) */
void finalise_spdist_matrix (spdist_matrix dist);
void complete_missing_spdist_from_global_spdist (spdist_matrix local, spdist_matrix global);
void copy_spdist_matrix_to_distance_matrix_upper (spdist_matrix spd, distance_matrix dist, bool use_means);
void del_spdist_matrix (spdist_matrix dist);


/*! \brief updates distances between species based on genes and gene-to-species mapping, with min on upper and mean on lower diagonal  */
void fill_species_dists_from_gene_dists (distance_matrix spdist, distance_matrix gendist, int *sp_id, bool use_upper_gene);
/*! \brief update global (over loci) species distances besed on local (within locus) species distances */
void update_species_dists_from_spdist (distance_matrix global, distance_matrix local, int *spexist);

int prepare_spdistmatrix_from_gene_species_map (spdist_matrix spdist, int *sp_id, int n_sp_id);
void fill_spdistmatrix_from_gene_dists (spdist_matrix spdist, distance_matrix gendist, int *sp_id, bool use_upper_gene);
/*! \brief initialise spdist_matrix with patristic distances from gdist vector of size n_gdist (1D) */
void fill_spdistmatrix_from_gene_dist_vector (spdist_matrix spdist, double *gdist, int n_gdist, int *sp_id);
void update_spdistmatrix_from_spdistmatrix (spdist_matrix global, spdist_matrix local);

#endif
