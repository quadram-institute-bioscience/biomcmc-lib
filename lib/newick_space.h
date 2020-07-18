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

/*! \file newick_space.h
 *  \brief Reads a list of trees in newick format and creates vector of topologies.
 *
 *  Currently does not check for duplicated trees, or same leaf names on a tree */

#ifndef _biomcmc_newick_space_h_
#define _biomcmc_newick_space_h_

#include "read_newick_trees.h" 

typedef struct newick_space_struct* newick_space;

/*! \brief Collection of topologies from tree file. Each topology will have its own char_vector */ 
struct newick_space_struct 
{
  int ntrees;      /*! \brief Number of trees originally in nexus file and compacted (only distinct topologies). */
  topology *t;     /*! \brief Vector of trees originally in nexus file and compacted. */
  int ref_counter; /*! \brief How many variables point to this structure (to see if memory can be free'd or not) */
};

newick_space new_newick_space (void);
void del_newick_space (newick_space nwk);
/*! \brief Convenience function to read one newick tree from file, skipping checks (comments, multiline trees, etc.) */
topology new_single_topology_from_newick_file (char *filename); 
newick_space new_newick_space_from_file (char *filename);
void update_newick_space_from_file (newick_space nwk, char *filename);
void update_newick_space_from_string (newick_space nwk, char *tree_string, size_t string_size);
void update_newick_space_from_topology (newick_space nwk, topology topol);

#endif
