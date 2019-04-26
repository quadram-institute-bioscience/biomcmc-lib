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

/*! \file topology_space.h
 *  \brief Reads tree files in nexus format and creates a vector of topologies. 
 *
 *  The #topology_space_struct is distinct from the #newick_space_struct since all trees here must share same
 *  char_vector (newick spaces can have general, uncomparable trees), and also since we store the distribution of trees
 *  (that is, each tree will have a frequency/representativity associated to it, as typical from Bayesian analyes).  */

#ifndef _biomcmc_topology_space_h_
#define _biomcmc_topology_space_h_

#include "newick_space.h" 

typedef struct topology_space_struct* topology_space;

/*! \brief Collection of topologies from tree file. When topologies have no branch lengths we store only unique
 * topologies */
struct topology_space_struct 
{
  int ntrees, ndistinct; /*! \brief Number of trees originally in nexus file and compacted (only distinct topologies). */
  topology *tree, *distinct; /*! \brief Vector of trees originally in nexus file and compacted. */
  double *freq;              /*! \brief frequency of each distinct topology (add up to one) */
  char_vector taxlabel;      /*! \brief Taxon names. */
  hashtable taxlabel_hash;   /*! \brief Lookup table with taxon names. */
  bool is_rooted;            /*! \brief If trees are unrooted, then branch lengths must be accounted for in some comparisons */ 
  char *filename;            /*! \brief name (without extension) of the originating file from where topology_space was read */
};

/*! \brief Read tree in newick format until char string_size, returning updated topolgy_space. Auxiliary for python module */
void add_string_with_size_to_topology_space (topology_space *tsp_address, char *long_string, size_t string_size, bool use_root_location);
/*! \brief Add topology to topology_space only if unrooted version is distinct, updating freqs, trees[] etc. Aux for python module */
void add_topology_to_topology_space_if_distinct (topology topol, topology_space tsp, double tree_weight, bool use_root_location);

/*! \brief Read tree file and store info in topology_space_struct with possible external hashtable to impose the leaf ordering. */
topology_space read_topology_space_from_file (char *seqfilename, hashtable external_taxhash, bool use_root_location);
/*! \brief lower level function where we can specify burnin and thinning factor, in iterations */ 
topology_space read_topology_space_from_file_with_burnin_thin (char *seqfilename, hashtable external_taxhash, int burnin, int thin, bool use_root_location);
/*!  \brief merge trees from two topology_space objects, assuming names hashtable is the same */ 
void merge_topology_spaces (topology_space ts1, topology_space ts2, double weight_ts1, bool use_root_location);
void  sort_topology_space_by_frequency(topology_space tsp, double *external_freqs) ; // INCOMPLETE

/*! \brief Save topology_space to a file, in format nexus w/ trprobs, up to "credible" cummul frequency */
void save_topology_space_to_trprobs_file (topology_space tsp, char *filename, double credible);
/*! \brief Quickly counts the number of leaves in a tree file, without storing any info. Assumes file and trees are well-formed */
int estimate_treesize_from_file (char *seqfilename);

/*! \brief Allocates memory for topology_space_struct (set of trees present in nexus file).  */
topology_space new_topology_space (void);
/*! \brief Free memory from topology_space_struct. */
void del_topology_space (topology_space tsp);

#endif
