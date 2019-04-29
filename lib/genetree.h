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

/*! \file genetree.h 
 *  \brief gene tree and species tree structures, for reconciliation etc. This is the high-level file with globally
 *  exposed functions/structures. 
 */

#ifndef _biomcmc_genetree_h_
#define _biomcmc_genetree_h_

#include "upgma.h"

typedef struct genetree_struct* genetree;
typedef struct speciestree_struct* speciestree;
typedef struct reconciliation_struct* reconciliation;
typedef struct splitset_struct* splitset;

struct genetree_struct
{
  topology t;
  reconciliation rec;
  speciestree sptre;  
  splitset split;
  int *distance, *minmax;
  int ref_counter;
};

struct speciestree_struct
{
  topology t;
  topol_node *mrca; /*! \brief triangular matrix of topol_nodes (LCA between topol_node::id (i-1) and j) in one dimension */
  int *spnames_order; /*! \brief Length+lexico order of sptree leaf names (not used unless added by user, when arbitrary leaf ordering is requested) */
  int ref_counter;
};

/*! \brief mapping between gene tree nodes (this) and (external) species tree nodes */
struct reconciliation_struct
{
  topol_node *map_d; /*! \brief Mapping of all nodes from gene to species (the first gene::nnodes are fixed) */
  topol_node *map_u; /*! \brief Mapping of all nodes from gene to species, assuming gene tree is upside down (unrooted TESTING version) */
  int *sp_id,    /*! \brief mapping of gene (this tolopogy) leaf to ID of taxon in species tree */
      *sp_count, /*! \brief how many copies of each species are present in this gene (used by deepcoal) */
      sp_size,   /*! \brief effective number of species present in gene family */
      size_diff, /*! \brief twice the difference in number of leaves between gene tree and (reduced/effective) species tree */
      *dup,      /*! \brief indexes of duplication nodes on gene tree, and number of such nodes (unused for now) */
      *ndup_d,   /*! \brief number of duplications below node (edge above node, since struct assume node == edge above it) */
      *ndup_u,   /*! \brief number of duplications above node (edge upside down, thus "children" are 'up' and 'sister') */
      *nlos_d,   /*! \brief number of losses below node and edge above node */  
      *nlos_u,   /*! \brief number of losses above node, including edge above it */
      ndups,     /*! \brief minimum number of duplications over all possible rootings, acc. to reconciliation_struct::dup */
      nloss,     /*! \brief number of losses corresponding to rooting (edge) that minimizes duplications */
      ndcos;     /*! \brief total number of deep coalescences (from nloss - 2 X ndups + size_diff) */ 
};

struct splitset_struct
{
  int size, spsize, spr, spr_extra, rf, hdist, hdist_reduced; /*! \brief spr, extra prunes for spr, rf distances and hdist=assignment cost */
  int n_g, n_s, n_agree, n_disagree;
  bipartition *g_split, *s_split, *agree, *disagree, *sp0; /* sp0 points to vec[0], s_split points to vec[x] */
  bipartition prune;
  hungarian h; /* hungarian method for solving the assignment between edges */
  bool match;  /*! \brief do we want to calculate the minimum cost assignment */
};

// convenience function
genetree new_genetree_speciestree_pair (topology gene, topology species);
/*! \brief Allocate space for new genetree_struct, given a gene topology and a specestree_struct */
genetree new_genetree (topology gene, speciestree sptre);
void del_genetree (genetree gtre);
/*! \brief Allocate space for new speciestree_struct, given a species topology and optionally the order of  species names */
speciestree new_speciestree (topology species, int *order_of_species_names);
void del_speciestree (speciestree sptre);
/*! \brief calculates all (discrete) distances and update min and max */
void genetree_speciestree_distances (genetree gtre, speciestree sptre);
/*! \brief from gene-species map index, count number of distinct species represented */
int count_species_in_index_species_gene (int *sp_id, int max_sp, int n_sp_id);

/*! \brief <debug function> dups.loss, ils calculation; accepts unseen speciestree_struct (i.e. updates mrca and pointers). Calls low-level hidden function. */
void genetree_reconcile_speciestree (genetree gtre, speciestree sptre);
/*! \brief <debug function> dSPR (level > 1), hdist (level > 0), and RF distances; doesn't need to update sptree pointer */
void genetree_dSPR_speciestree (genetree gtre, speciestree sptre, int level);

#endif
