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

/*! \file read_newick_trees.h
 *  \brief Low-level functions for reading newick strings 
 *
 *  Currently does not check for duplicated trees, or same leaf names on a tree */

#ifndef _biomcmc_read_newick_trees_h
#define _biomcmc_read_newick_trees_h

#include "topology_common.h" 
#include "nexus_common.h" 

typedef struct newick_node_struct* newick_node;
typedef struct newick_tree_struct* newick_tree; /*! \brief newick trees have minimal information, unlike topology_struct */

struct newick_node_struct
{
  newick_node up, right, left; /*! \brief Parent and children nodes. */
  int id;               /*! \brief Initial pre-order numbering of node. */
  double branch_length; /*! \brief Branch length from node to node->up. */
  char *taxlabel;       /*! \brief Leaf sequence name */
};

struct newick_tree_struct
{
  newick_node *nodelist; /*! \brief Vector with pointers to every internal node. */
  newick_node *leaflist; /*! \brief Vector with pointers to tree leaves. */
  newick_node root;      /*! \brief Pointer to root node. */
  bool has_branches;     /*! \brief Boolean saying if tree has branch lengths or not. (topology alsways has, even if one-zero) */
  int nnodes, nleaves;   /*! \brief Number of nodes (including leaves), and number of leaves */
};

/*! \brief Allocates memory for newick_tree_struct. */
newick_tree new_newick_tree (int nleaves);
/*! \brief Frees memory used by tree. */
void del_newick_tree (newick_tree T);
/*! \brief Copy information from newick_tree struct to topology_struct; newick_space copies taxlabels but topology_space
 * (from nexus files) share the taxlabel and thus don't copy from #newick_tree_struct  */
void copy_topology_from_newick_tree (topology tree, newick_tree nwk_tree, bool create_tree_taxlabel);
/*! \brief Creates newick_tree structure. */ 
newick_tree new_newick_tree_from_string (char *external_string);
/*! \brief Recursive function that creates a node based on parenthetic structure. */
newick_node subtree_newick_tree (newick_tree tree, char *lsptr, char *rsptr, int *node_id, newick_node up);
/*! \brief Counts the number of leaves and resolves (one) trifurcation of tree string. */
int number_of_leaves_in_newick (char **string);

#endif
