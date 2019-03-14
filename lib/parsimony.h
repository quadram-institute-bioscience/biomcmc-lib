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

/*! \file hashtable.h
 *  \brief double hashing open-address hash table using strings as key -- also has distance matrix, that can be used in
 *  alignments and trees
 *
 * Hash tables allow us to search for the position of a key (taxa name) without scanning the whole vector (like in 
 * sequencial search). This code is derived from the software DCM3, released under the GPL license (Copyright (C) 
 * 2004 The University of Texas at Austin). 
 */

#ifndef _biomcmc_parsimony_mrp_
#define _biomcmc_parsimony_mrp_

#include "topology_common.h"

typedef struct binary_mrp_matrix_struct* binary_mrp_matrix; 
typedef struct mrp_parsimony_struct* mrp_parsimony;

/*! \brief matrix representation with parsimony (01 10 11 sequences) */
struct binary_mrp_matrix_struct {
  int ntax, nchar, npat;  /*!< \brief number of taxa, sites, and distinct patterns */
  bool **s;               /*!< \brief 1 (01) and 2 (10) are the two binary states, with 3 (11) being undetermined */
  int *pattern_freq;      /*!< \brief frequency of pattern. */
  int ref_counter;        /*!< \brief how many places have a pointer to this instance */
};

struct mrp_parsimony_struct {
  int *score;      /*!< \brief parsimony score per pattern */
  binary_mrp_matrix external, internal; /*!< \brief binary matrices for leaves and for internal nodes */
  int ref_counter; /*!< \brief how many places have a pointer to this instance */
};

binary_mrp_matrix new_binary_mrp_matrix (int n_sequences, int n_sites);
void del_binary_mrp_matrix (binary_mrp_matrix mrp);
mrp_parsimony new_mrp_parsimony (int n_sequences, int n_sites);
void del_mrp_parsimony (mrp_parsimony pars);
/*! \brief given a map[] with location in sptree of gene tree leaves, update binary matrix with splits from genetree */
void update_binary_mrp_matrix_from_topology (binary_mrp_matrix mrp, topology t, int *map);

#endif
