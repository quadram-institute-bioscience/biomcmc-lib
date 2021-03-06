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

/*! \file reconciliation.h 
 *  \brief low-level file for gene tree and species tree reconciliation. This file is hidden from user and contains the
 *  LCA-based reconciliation distances.  */

#ifndef _biomcmc_reconciliation_h_
#define _biomcmc_reconciliation_h_

#include "genetree.h"  // genetree.c -> genetree.h, aux.h -> genetree.h 

/*! \brief Allocate space for new reconciliation_struct (other functions defined in topology_mrca.c) */
reconciliation new_reconciliation (int gene_nleaves, int sp_nleaves);
/*! \brief Create new reconciliation struct and copy values (except species tree info) from another struct */
reconciliation new_reconciliation_from_reconciliation (int gene_nleaves, int sp_nleaves, reconciliation from);
/*! \brief release allocated memory for reconciliation_struct */
void del_reconciliation (reconciliation r);

/*! \brief Fill rec->sp_count[] with the number of representatives of each species (idexed by rec->sp_id[]) */
void initialize_reconciliation_sp_count (reconciliation rec, int n_sp, int n_idx);

/*! \brief transform indexes found in index_sptaxa_to_genetaxa() to pointers to species nodes */
void initialize_reconciliation_from_new_species_tree (genetree gtre, speciestree sptre);

/*! \brief Find reconciliation map between gene and species trees */
void reconciliation_gene_tree_reconcile (genetree gtre, speciestree sptre);

#endif
